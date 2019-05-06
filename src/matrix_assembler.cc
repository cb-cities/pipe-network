
#include <cmath>
#include <matrix_assembler.h>

// Constructor
pipenetwork::MatrixAssembler::MatrixAssembler(bool pdd_mode) {
  pdd_ = pdd_mode;
  jac_ = std::make_shared<Eigen::SparseMatrix<double>>();
  variable_vec_ = std::make_shared<Eigen::VectorXd>();
  residual_vec_ = std::make_shared<Eigen::VectorXd>();
};

// Obtain global nodal and pipe indices and pointers from meshes
void pipenetwork::MatrixAssembler::global_nodal_pipe_indices(
    const std::shared_ptr<pipenetwork::Mesh>& mesh) {
  global_nodes_.clear();
  global_pipes_.clear();
  for (const auto& node : mesh->nodes_) global_nodes_.emplace(node);
  for (const auto& pipe : mesh->pipes_) global_pipes_.emplace(pipe);
  nnode_ = global_nodes_.size();
  npipe_ = global_pipes_.size();
}

// Initialize variable vector
// 1 to nnode element: nodal head
// nnode+1 to 2*nnode element: nodal demand
// 2*nnode+1 to 2*nnode+npipe element: pipe discharge
void pipenetwork::MatrixAssembler::assemble_variable_vector() {

  variable_vec_->resize(2 * nnode_ + npipe_);

  for (const auto& node : global_nodes_) {
    Index index_nh = node.first;
    Index index_nd = node.first + nnode_;

    // nodal head vector
    variable_vec_->coeffRef(index_nh) = node.second->head();

    // nodal demand vector
    variable_vec_->coeffRef(index_nd) = node.second->demand();
  }

  // pipe discharge vector
  for (const auto& pipe : global_pipes_) {
    Index index_pd = pipe.first + 2 * nnode_;
    variable_vec_->coeffRef(index_pd) = pipe.second->discharge();
  }
}

// assember the pressure demand part of jacobian matrix
std::vector<Eigen::Triplet<double>>
    pipenetwork::MatrixAssembler::construct_demand_jac(
        const std::shared_ptr<pipenetwork::Node>& node, Index index) {
  std::vector<Eigen::Triplet<double>> update;
  if (not node->isres()) {
    // construct jacH part
    update.emplace_back(nnode_ + npipe_ + index, nnode_ + index, 1);
    // construct jacG part
    if (pdd_) {
      update.emplace_back(nnode_ + npipe_ + index, index,
                          get_pressure_head_jacob(node));
    }
  } else {
    // JacH will be 0, consturct jac G
    update.emplace_back(nnode_ + npipe_ + index, index, 1);
  }
  return update;
}

// Assemble Jacobian matrix
//                 nodal_head  nodal_discharge(demand) pipe_discharge
// nodal_balance   sub_jac_A        sub_jac_B         sub_jac_C
// headloss        sub_jac_D        sub_jac_E         sub_jac_F
// demand_pressure sub_jac_G        sub_jac_H         sub_jac_I
//
// Nodal balance equation:
// Residual = flow = Inflow - Outflow - Demand
//
// Headloss equation (Hazen-Williams):
// Residual = (Start_node_head - end_node_head) - \frac{10.67 \times length
//            \times pipe_discharge^1.852}{pipe_roughness^1.852 \times
//            (2radius)^4.8704}
//
// demand_pressure equation :
// Residual = demand-demand (demand driven)
// Residual = pressure_demand_func-demand (demand driven)
//
// sub_jac_A: Derivative of nodal balance equation with respect to nodal head,
//            0.
// sub_jac_B: Derivative of nodal balance equation with respect to nodal
//            discharge (demand), -1.
// sub_jac_C: Derivative of nodal balance equation with respect to pipe
//            discharge, 1 for inlink, -1 for out link.
// sub_jac_D: Derivative of headloss equation with respect to nodal head, 1 for
//            start node, -1 for end node.
// sub_jac_E: Derivative of headloss equation with respect to nodal discharge,
//            0.
// sub_jac_F: Derivative of headloss equation with respect to pipe discharge,
//            for corresponding pipe, \frac{-1.852 \times 10.67 \times
//            length}{pow(pipe_roughness,1.852) \times pow(2radius,4.8704)}
//            \times pow(pipe_discharge,0.852).
// sub_jac_G: Derivative of demand_pressure equation with respect to nodal head,
//            0 for junctions, 1 for reservoir/tank of demand driven model,
//            nonlinear for pressure_demand model.
// sub_jac_H: Derivative of demand_pressure equation with respect to nodal
// discharge (demand),
//            1 for junctions, 0 for reservoir/tank.
// sub_jac_I: Derivative of demand_pressure equation with respect to pipe
// discharge,
//            0.
//
// Each node has one nodal balance equation, each pipe has one headloss
// equation, Thus the Jacobian has (2*nnode+npipe) row and (2*nnode+npipe)
// column
//

void pipenetwork::MatrixAssembler::assemble_jacobian() {
  // Check network
  if (nnode_ <= 0 || npipe_ <= 0)
    throw std::runtime_error(
        "No node or pipe index pairs created to assemble Jacobian matrix");

  std::vector<Eigen::Triplet<double>> update;
  update.reserve(nnode_ + 5 * npipe_);

  // Iterate through all nodes
  for (const auto& node : global_nodes_) {
    Index index = node.first;
    // construct jacB part
    update.emplace_back(index, nnode_ + index, -1);
    // construct jacH, jacG part
    auto demand_update = construct_demand_jac(node.second, index);
    update.insert(std::end(update), std::begin(demand_update),
                  std::end(demand_update));
  }

  // Iterate through all pipes
  for (const auto& pipe : global_pipes_) {
    Index index_pipe = pipe.first;

    Index index_node1 = pipe.second->nodes().at(0)->id();
    Index index_node2 = pipe.second->nodes().at(1)->id();

    // construct jacC part
    if (pipe.second->iter_discharge() >= 0) {
      update.emplace_back(index_node1, 2 * nnode_ + index_pipe, -1);
      update.emplace_back(index_node2, 2 * nnode_ + index_pipe, 1);
    } else {
      update.emplace_back(index_node1, 2 * nnode_ + index_pipe, 1);
      update.emplace_back(index_node2, 2 * nnode_ + index_pipe, -1);
    }
    // construct jacD part
    update.emplace_back(nnode_ + index_pipe, index_node1, 1);
    update.emplace_back(nnode_ + index_pipe, index_node2, -1);
    // construct jacF part (Hazen-Williams)
    double deriv_pipe_discharge = pipe.second->deriv_hazen_williams_discharge();
    update.emplace_back(nnode_ + index_pipe, 2 * nnode_ + index_pipe,
                        deriv_pipe_discharge);
  }

  jac_->resize(2 * nnode_ + npipe_, 2 * nnode_ + npipe_);
  jac_->setFromTriplets(update.begin(), update.end());
}

// Calculate and assemble residual vector
// Nodal balance equation:
// Residual = flow = Inflow - Outflow - Demand
//
//
// Headloss equation (Hazen-Williams):
// Residual = (Start_node_head - end_node_head) - \frac{10.67 \times length
//            \times pipe_discharge^1.852}{pipe_roughness^1.852 \times
//            (2radius)^4.8704}

// demand_pressure equation :
// Residual = demand_var-demand_desired (demand driven mode)
// Residual = demand_var - pressure_demand_func(pressure demand driven mode)
// The network has nnode Nodal balance equation and npipe Headloss equation
// Thus the residual vector has (nnode+npipe) elements
void pipenetwork::MatrixAssembler::assemble_residual_vector() {
  residual_vec_->resize(2 * nnode_ + npipe_);
  residual_vec_->setZero();

  // Iterate through all pipes
  for (auto& pipe : global_pipes_) {
    Index index_pipe = pipe.first;
    Index index_node1 = pipe.second->nodes().at(0)->id();
    Index index_node2 = pipe.second->nodes().at(1)->id();
    // Calculate headloss residual
    pipe.second->compute_headloss_hazen_williams();
    double head1 = 0;
    double head2 = 0;
    head1 = pipe.second->nodes().at(0)->head();
    head2 = pipe.second->nodes().at(1)->head();

    auto head_res = (-1) * ((head1 - head2) - pipe.second->headloss());
    residual_vec_->coeffRef(nnode_ + index_pipe) = head_res;

    // Calculate the nodal balance residual (in and out flow part)
    residual_vec_->coeffRef(index_node1) +=
        (-1) * (-1) * pipe.second->iter_discharge();
    residual_vec_->coeffRef(index_node2) +=
        (-1) * pipe.second->iter_discharge();
  }

  // Iterate through all nodes
  for (const auto& node : global_nodes_) {
    Index index = node.first;
    // Calculate the nodal balance residual (demand part)
    if (node.second->isdischarge())
      residual_vec_->coeffRef(index) -= (-1) * node.second->iter_demand();
    // Calculate the demand pressure residual
    // Case 1 the node is Reservior/tank
    if (node.second->isres()) {
      residual_vec_->coeffRef(nnode_ + npipe_ + index) =
          node.second->head() - node.second->elevation();
    }
    // case 2 the node is a junction
    else {
      if (not pdd_) {
        residual_vec_->coeffRef(nnode_ + npipe_ + index) =
            node.second->iter_demand() - node.second->demand();
      } else {
        assemble_pdd_residual(node.second, index);
      }
    }
  }
}

// Apply variables (head and discharge) back to nodes and pipes
void pipenetwork::MatrixAssembler::apply_variables() {

  // Iterate through nodes, assign nodal variables
  for (auto& node : global_nodes_) {
    // Get index of target variables
    Index index_nh = node.first;
    // Assign head
    node.second->head(variable_vec_->coeff(index_nh));
    // Assign demand
    node.second->iter_demand(variable_vec_->coeff(nnode_ + index_nh));
  }

  // Iterate through pipes, assign pipe discharge during iteration
  for (auto& pipe : global_pipes_) {
    // Get index of target variables
    Index index_pd = pipe.first + 2 * nnode_;
    // Assign discharge
    pipe.second->iter_discharge(variable_vec_->coeff(index_pd));
  }
}

void pipenetwork::MatrixAssembler::assemble_pdd_residual(
    const std::shared_ptr<pipenetwork::Node>& node, Index index) {
  auto pressure = node->head() - node->elevation();
  double res;

  // case 1: pressure below min pressure, no discharge
  if (pressure < node->min_pressure()) {
    res = -1 * (node->iter_demand() - node->demand() * node->pdd_slope() *
                                          (pressure - node->min_pressure()));
  }
  // case 2: pressure just above min pressure, use polynomial to approximate
  // nonlinear pressure-demand function for stabilization considerations
  else if ((pressure > node->min_pressure()) and
           (pressure <= (node->min_pressure() + node->pdd_smooth_delta()))) {
    auto pdd_coeff_1 = node->get_pdd_poly_coef_1();
    res = -1 * (node->iter_demand() -
                node->demand() *
                    (pdd_coeff_1[0] * std::pow(pressure, 3) +
                     pdd_coeff_1[1] * std::pow(pressure, 2) +
                     pdd_coeff_1[2] * std::pow(pressure, 1) + pdd_coeff_1[3]));
  }
  // case 3: pressure-demand nonlinear function
  else if ((pressure > (node->min_pressure() + node->pdd_smooth_delta())) and
           (pressure <= (node->norm_pressure() - node->pdd_smooth_delta()))) {
    res = -1 * (node->iter_demand() -
                node->demand() *
                    std::pow((pressure - node->min_pressure()) /
                                 (node->norm_pressure() - node->min_pressure()),
                             0.5));
  }
  // case 4:pressure close to normal pressure, use polynomial to approximate
  // nonlinear pressure-demand function for stabilization considerations
  else if ((pressure > ((node->norm_pressure() - node->pdd_smooth_delta())) and
            (pressure <= node->norm_pressure()))) {
    auto pdd_coeff_2 = node->get_pdd_poly_coef_2();
    res = -1 * (node->iter_demand() -
                node->demand() *
                    (pdd_coeff_2[0] * std::pow(pressure, 3) +
                     pdd_coeff_2[1] * std::pow(pressure, 2) +
                     pdd_coeff_2[2] * std::pow(pressure, 1) + pdd_coeff_2[3]));
  }
  // case 5: pressure above normal pressure, demand can be met
  else {
    res = -1 * (node->iter_demand() - node->demand());
  }
  // assign the residual
  residual_vec_->coeffRef(nnode_ + npipe_ + index) = res;
}

double pipenetwork::MatrixAssembler::get_pressure_head_jacob(
    const std::shared_ptr<pipenetwork::Node>& node) {

  auto pressure = node->head() - node->elevation();
  double res;

  // case 1: pressure below min pressure, no discharge
  if (pressure < node->min_pressure()) {
    res = 0;
  }
  // case 2: pressure just above min pressure, use polynomial to approximate
  // nonlinear pressure-demand function for stabilization considerations
  else if ((pressure > node->min_pressure()) and
           (pressure <= (node->min_pressure() + node->pdd_smooth_delta()))) {
    auto pdd_coeff_1 = node->get_pdd_poly_coef_1();
    res = -1 * node->demand() *
          (3 * pdd_coeff_1[0] * std::pow(pressure, 2) +
           2 * pdd_coeff_1[1] * std::pow(pressure, 1) + pdd_coeff_1[2]);
  }
  // case 3: pressure-demand nonlinear function
  else if ((pressure > (node->min_pressure() + node->pdd_smooth_delta())) and
           (pressure <= (node->norm_pressure() - node->pdd_smooth_delta()))) {
    res = -1 * (0.5 * node->demand() *
                std::pow((pressure - node->min_pressure()) /
                             (node->norm_pressure() - node->min_pressure()),
                         -0.5));
  }
  // case 4:pressure close to normal pressure, use polynomial to approximate
  // nonlinear pressure-demand function for stabilization considerations
  else if ((pressure > ((node->norm_pressure() - node->pdd_smooth_delta())) and
            (pressure <= node->norm_pressure()))) {
    auto pdd_coeff_2 = node->get_pdd_poly_coef_2();
    res = -1 * node->demand() *
          (3 * pdd_coeff_2[0] * std::pow(pressure, 2) +
           2 * pdd_coeff_2[1] * std::pow(pressure, 1) + pdd_coeff_2[2]);
  }
  // case 5: pressure above normal pressure, demand can be met
  else {
    res = 0;
  }

  return res;
}
