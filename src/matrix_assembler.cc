
#include <matrix_assembler.h>

#include "matrix_assembler.h"

// Constructor
pipenetwork::MatrixAssembler::MatrixAssembler() {
  jac_ = std::make_shared<Eigen::SparseMatrix<double>>();
  variable_vec_ = std::make_shared<Eigen::VectorXd>();
  residual_vec_ = std::make_shared<Eigen::VectorXd>();
  id_map_ = std::make_shared<std::map<Index, Index>>();
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

// Obtain global nodal and pipe indices and pointers from meshes
void pipenetwork::MatrixAssembler::global_nodal_pipe_indices(
    const std::shared_ptr<pipenetwork::Mesh>& mesh,
    const unsigned nnode_known) {
  nnode_known_ = nnode_known;
  global_nodal_pipe_indices(mesh);
}

// Initialize variable vector
// 1 to nnode element: nodal head
// nnode+1 to 2*nnode element: nodal discharge
// 2*nnode+1 to 2*nnode+npipe element: pipe discharge
void pipenetwork::MatrixAssembler::assemble_variable_vector() {

  variable_vec_->resize(2 * nnode_ + npipe_);

  for (const auto& node : global_nodes_) {
    Index index_nh = node.first;
    Index index_nd = node.first + nnode_;

    // nodal head vector
    // If head of the ndoe is unknown (hasn't been assigned),
    // initialize to zero
    if ((node.second)->ishead())
      variable_vec_->coeffRef(index_nh) = node.second->head();
    else
      variable_vec_->coeffRef(index_nh) = 0.0;

    // nodal discharge vector
    // If discharge of the node is unknown (hasn't been assigned),
    // initialize to zero
    if ((node.second)->isdischarge())
      variable_vec_->coeffRef(index_nd) = node.second->discharge();
    else
      variable_vec_->coeffRef(index_nd) = 0.0;
  }

  // pipe discharge vector
  // Calculated according to nodal heads at two end and Hazen-Williams
  // equation. If any of the nodal head is unknown, initialize the
  // discharge to 0.001
  for (const auto& pipe : global_pipes_) {
    Index index_pd = pipe.first + 2 * nnode_;
    variable_vec_->coeffRef(index_pd) = pipe.second->discharge();
  }
}

// Apply variables (head and discharge) back to nodes and pipes
void pipenetwork::MatrixAssembler::apply_variables() {

  // Iterate through nodes, assign nodal variables
  for (auto& node : global_nodes_) {
    // Get index of target variables
    Index index_nh = node.first;
    Index index_nd = node.first + nnode_;
    // Assign head
    node.second->head(variable_vec_->coeff(index_nh));
    // Assign discharge
    node.second->discharge(variable_vec_->coeff(index_nd));
  }

  // Iterate through pipes, assign pipe discharge during iteration
  for (auto& pipe : global_pipes_) {
    // Get index of target variables
    Index index_pd = pipe.first + 2 * nnode_;
    // Assign discharge
    pipe.second->iter_discharge(variable_vec_->coeff(index_pd));
  }
}

// Calculate and assemble residual vector
// Nodal balance equation:
// Residual = flow = Inflow - Outflow - Demand
// Headloss equation (Hazen-Williams):
// Residual = (Start_node_head - end_node_head) - \frac{10.67 \times length
//            \times pipe_discharge^1.852}{pipe_roughness^1.852 \times
//            (2radius)^4.8704}
// The network has nnode Nodal balance equation and npipe Headloss equation
// Thus the residual vector has (nnode+npipe) elements
void pipenetwork::MatrixAssembler::assemble_residual_vector() {
  residual_vec_->resize(nnode_ + npipe_);
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
    if (pipe.second->nodes().at(0)->ishead())
      head1 = pipe.second->nodes().at(0)->head();
    if (pipe.second->nodes().at(1)->ishead())
      head2 = pipe.second->nodes().at(1)->head();
    residual_vec_->coeffRef(nnode_ + index_pipe) =
        (-1) * ((head1 - head2) - pipe.second->headloss());
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
      residual_vec_->coeffRef(index) -= (-1) * node.second->discharge();
  }
}

// Assemble Jacobian matrix
//                 nodal_head    nodal_discharge    pipe_discharge
// nodal_balance   sub_jac_A        sub_jac_B         sub_jac_C
// headloss        sub_jac_D        sub_jac_E         sub_jac_F
//
// Nodal balance equation:
// Residual = flow = Inflow - Outflow - Demand
// Headloss equation (Hazen-Williams):
// Residual = (Start_node_head - end_node_head) - \frac{10.67 \times length
//            \times pipe_discharge^1.852}{pipe_roughness^1.852 \times
//            (2radius)^4.8704}
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
//
// Each node has one nodal balance equation, each pipe has one headloss
// equation, Thus the Jacobian has (nnode+npipe) row and (2*nnode+npipe) column
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

  jac_->resize(nnode_ + npipe_, 2 * nnode_ + npipe_);
  jac_->setFromTriplets(update.begin(), update.end());
}

//**************************************************************
//***** Make all nodal discharge known value, test purpose *****
//**************************************************************

// Initialize variable vector
// 1 to nnode element: nodal head
// nnode+1 to nnode+npipe element: pipe discharge
void pipenetwork::MatrixAssembler::sim_assemble_variable_vector() {

  variable_vec_->resize(nnode_ + npipe_);

  for (const auto& node : global_nodes_) {
    Index index_nh = node.first;

    // nodal head vector
    // If head of the ndoe is unknown (hasn't been assigned),
    // initialize to zero
    if ((node.second)->ishead())
      variable_vec_->coeffRef(index_nh) = node.second->head();
    else
      variable_vec_->coeffRef(index_nh) = 0.0;
  }

  // pipe discharge vector
  // Calculated according to nodal heads at two end and Hazen-Williams
  // equation. If any of the nodal head is unknown, initialize the
  // discharge to 0.001
  for (const auto& pipe : global_pipes_) {
    Index index_pd = pipe.first + nnode_;
    variable_vec_->coeffRef(index_pd) = pipe.second->iter_discharge();
  }
}

// Apply variables (head and discharge) back to nodes and pipes
void pipenetwork::MatrixAssembler::sim_apply_variables() {

  // Iterate through nodes, assign nodal variables
  for (auto& node : global_nodes_) {
    // Get index of target variables
    Index index_nh = node.first;
    // Assign head
    node.second->head(variable_vec_->coeff(index_nh));
  }

  // Iterate through pipes, assign pipe discharge during iteration
  for (auto& pipe : global_pipes_) {
    // Get index of target variables
    Index index_pd = pipe.first + nnode_;
    // Assign discharge
    pipe.second->iter_discharge(variable_vec_->coeff(index_pd));
  }
}

// Assemble Jacobian matrix
//                 nodal_head    pipe_discharge
// nodal_balance   sub_jac_A        sub_jac_C
// headloss        sub_jac_D        sub_jac_F
//
// Nodal balance equation:
// Residual = flow = Inflow - Outflow - Demand
// Headloss equation (Hazen-Williams):
// Residual = (Start_node_head - end_node_head) - \frac{10.67 \times length
//            \times pipe_discharge^1.852}{pipe_roughness^1.852 \times
//            (2radius)^4.8704}
//
// sub_jac_A: Derivative of nodal balance equation with respect to nodal head,
//            0.
// sub_jac_C: Derivative of nodal balance equation with respect to pipe
//            discharge, 1 for inlink, -1 for out link.
// sub_jac_D: Derivative of headloss equation with respect to nodal head, 1 for
//            start node, -1 for end node.
// sub_jac_F: Derivative of headloss equation with respect to pipe discharge,
//            for corresponding pipe, \frac{-1.852 \times 10.67 \times
//            length}{pow(pipe_roughness,1.852) \times pow(2radius,4.8704)}
//            \times pow(pipe_discharge,0.852).
//
// Each node has one nodal balance equation, each pipe has one headloss
// equation, Thus the Jacobian has (nnode+npipe) row and (nnode+npipe) column
//
void pipenetwork::MatrixAssembler::sim_assemble_jacobian() {
  // Check network
  if (nnode_ <= 0 || npipe_ <= 0)
    throw std::runtime_error(
        "No node or pipe index pairs created to assemble Jacobian matrix");

  std::vector<Eigen::Triplet<double>> update;
  update.reserve(5 * npipe_);

  // Iterate through all pipes
  for (const auto& pipe : global_pipes_) {
    Index index_pipe = pipe.first;

    Index index_node1 = pipe.second->nodes().at(0)->id();
    Index index_node2 = pipe.second->nodes().at(1)->id();

    // construct jacC part
    if (pipe.second->iter_discharge() >= 0) {
      update.emplace_back(index_node1, nnode_ + index_pipe, -1);
      update.emplace_back(index_node2, nnode_ + index_pipe, 1);
    } else {
      update.emplace_back(index_node1, nnode_ + index_pipe, 1);
      update.emplace_back(index_node2, nnode_ + index_pipe, -1);
    }
    // construct jacD part
    update.emplace_back(nnode_ + index_pipe, index_node1, 1);
    update.emplace_back(nnode_ + index_pipe, index_node2, -1);
    // construct jacF part (Hazen-Williams)
    double deriv_pipe_discharge = pipe.second->deriv_hazen_williams_discharge();
    update.emplace_back(nnode_ + index_pipe, nnode_ + index_pipe,
                        deriv_pipe_discharge);
  }

  jac_->resize(nnode_ + npipe_, nnode_ + npipe_);
  jac_->setFromTriplets(update.begin(), update.end());
}

// Initialize variable vector
// 1 to nnode-nnode_known element: nodal head
// nnode-nnode_known+1 to nnode-nnode_known+npipe element: pipe discharge
void pipenetwork::MatrixAssembler::sim_assemble_variable_vector_v2() {
  variable_vec_->resize(nnode_ + npipe_ - nnode_known_);

  Index index_var = 0;
  for (const auto& node : global_nodes_) {
    Index index_nh = node.first;

    // nodal head vector
    // If head of the node is known (has been assigned),
    // skip it to avoid overconstrain
    // else initialize it to be zero
    if ((node.second)->ishead()) {
      continue;
    } else {
      variable_vec_->coeffRef(index_var) = 0.0;
      id_map_->emplace(index_nh, index_var);
    }
    index_var += 1;
  }

  // pipe discharge vector
  // Get the initialized pipe discharge from the mesh
  for (const auto& pipe : global_pipes_) {
    Index index_pd = pipe.first + nnode_ - nnode_known_;
    variable_vec_->coeffRef(index_pd) = pipe.second->discharge();
  }
}

// Apply variables (head and discharge) back to nodes and pipes
void pipenetwork::MatrixAssembler::sim_apply_variables_v2() {
  // Iterate through nodes, assign nodal variables
  for (auto& node : global_nodes_) {
    // Get index of target variables
    Index index_nh = node.first;

    // Assign variables back to mesh
    if (check_node_avail(index_nh)) {
      // found
      Index index_var = id_map_->at(index_nh);
      // Assign head
      node.second->head(variable_vec_->coeff(index_var));
    }
  }

  // Iterate through pipes, assign pipe discharge during iteration
  for (auto& pipe : global_pipes_) {
    // Get index of target variables
    Index index_pd = pipe.first + nnode_ - nnode_known_;
    // Assign discharge
    pipe.second->iter_discharge(variable_vec_->coeff(index_pd));
  }
}

// Check if the variable vector contains given node
bool pipenetwork::MatrixAssembler::check_node_avail(Index node_id) {
  if (id_map_->find(node_id) != id_map_->end()) {
    return true;
  }
  return false;
}

// Assemble Jacobian matrix (for unknown variables only)
//                 nodal_head    pipe_discharge
// nodal_balance   sub_jac_A        sub_jac_C
// headloss        sub_jac_D        sub_jac_F
//
// Nodal balance equation:
// Residual = flow = Inflow - Outflow - Demand
// Headloss equation (Hazen-Williams):
// Residual = (Start_node_head - end_node_head) - \frac{10.67 \times length
//            \times pipe_discharge^1.852}{pipe_roughness^1.852 \times
//            (2radius)^4.8704}
//
// sub_jac_A: Derivative of nodal balance equation with respect to nodal head,
//            0.
// sub_jac_C: Derivative of nodal balance equation with respect to pipe
//            discharge, 1 for inlink, -1 for out link.
// sub_jac_D: Derivative of headloss equation with respect to nodal head, 1 for
//            start node, -1 for end node.
// sub_jac_F: Derivative of headloss equation with respect to pipe discharge,
//            for corresponding pipe, \frac{-1.852 \times 10.67 \times
//            length}{pow(pipe_roughness,1.852) \times pow(2radius,4.8704)}
//            \times pow(pipe_discharge,0.852).
//
// Each node has one nodal balance equation, each pipe has one headloss
// equation, Thus the Jacobian has (nnode+npipe) row and (nnode+npipe) column
//
void pipenetwork::MatrixAssembler::sim_assemble_jacobian_v2() {
  // Check network
  if (nnode_ <= 0 || npipe_ <= 0)
    throw std::runtime_error(
        "No node or pipe index pairs created to assemble Jacobian matrix");

  std::vector<Eigen::Triplet<double>> update;
  update.reserve(5 * npipe_);

  // Iterate through all pipes
  for (const auto& pipe : global_pipes_) {
    Index index_pipe = pipe.first;

    Index index_node1 = pipe.second->nodes().at(0)->id();
    Index index_node2 = pipe.second->nodes().at(1)->id();

    // construct jacC  part
    if (pipe.second->iter_discharge() >= 0) {

      if (check_node_avail(index_node1)) {
        Index index_var = id_map_->at(index_node1);
        update.emplace_back(index_var, nnode_ - nnode_known_ + index_pipe, -1);
      }
      if (check_node_avail(index_node2)) {
        Index index_var = id_map_->at(index_node2);
        update.emplace_back(index_var, nnode_ - nnode_known_ + index_pipe, 1);
      }
    } else {
      if (check_node_avail(index_node1)) {
        Index index_var = id_map_->at(index_node1);
        update.emplace_back(index_var, nnode_ - nnode_known_ + index_pipe, 1);
      }
      if (check_node_avail(index_node2)) {
        Index index_var = id_map_->at(index_node2);
        update.emplace_back(index_var, nnode_ - nnode_known_ + index_pipe, -1);
      }
    }

    // construct jacD part
    if (check_node_avail(index_node1)) {
      Index index_var = id_map_->at(index_node1);
      update.emplace_back(nnode_ - nnode_known_ + index_pipe, index_var, 1);
    }
    if (check_node_avail(index_node2)) {
      Index index_var = id_map_->at(index_node2);
      update.emplace_back(nnode_ - nnode_known_ + index_pipe, index_var, -1);
    }

    // construct jacF part (Hazen-Williams)
    double deriv_pipe_discharge = pipe.second->deriv_hazen_williams_discharge();
    update.emplace_back(nnode_ - nnode_known_ + index_pipe,
                        nnode_ - nnode_known_ + index_pipe,
                        deriv_pipe_discharge);
  }

  jac_->resize(nnode_ + npipe_ - nnode_known_, nnode_ + npipe_ - nnode_known_);
  jac_->setFromTriplets(update.begin(), update.end());
}

// Calculate and assemble residual vector for unknown variables only
// Nodal balance equation:
// Residual = flow = Inflow - Outflow - Demand
// Headloss equation (Hazen-Williams):
// Residual = (Start_node_head - end_node_head) - \frac{10.67 \times length
//            \times pipe_discharge^1.852}{pipe_roughness^1.852 \times
//            (2radius)^4.8704}
// The network has nnode Nodal balance equation and npipe Headloss equation
// Thus the residual vector has (nnode+npipe) elements
void pipenetwork::MatrixAssembler::assemble_residual_vector_v2() {
  residual_vec_->resize(nnode_ + npipe_ - nnode_known_);
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
    if (pipe.second->nodes().at(0)->ishead())
      head1 = pipe.second->nodes().at(0)->head();
    if (pipe.second->nodes().at(1)->ishead())
      head2 = pipe.second->nodes().at(1)->head();
    residual_vec_->coeffRef(nnode_ - nnode_known_ + index_pipe) =
        (-1) * ((head1 - head2) - pipe.second->headloss());
    // Calculate the nodal balance residual (in and out flow part)
    if (check_node_avail(index_node1)) {
      Index index_var = id_map_->at(index_node1);
      residual_vec_->coeffRef(index_var) +=
          (-1) * (-1) * pipe.second->iter_discharge();
    }

    if (check_node_avail(index_node2)) {
      Index index_var = id_map_->at(index_node2);
      residual_vec_->coeffRef(index_var) +=
          (-1) * pipe.second->iter_discharge();
    }
  }

  // Iterate through all nodes
  for (const auto& node : global_nodes_) {
    Index index = node.first;
    // Calculate the nodal balance residual (demand part)
    if (check_node_avail(index)) {
      Index index_var = id_map_->at(index);
      if (node.second->isdischarge())
        residual_vec_->coeffRef(index_var) -= (-1) * node.second->discharge();
    }
  }
}