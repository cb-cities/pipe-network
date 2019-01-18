// Constructor
MatrixAssembler::MatrixAssembler() {
  jac_ = std::make_shared<Eigen::SparseMatrix<double>>();
  node_head_vec_ = std::make_shared<Eigen::VectorXd>();
  node_discharge_vec_ = std::make_shared<Eigen::VectorXd>();
  pipe_discharge_vec_ = std::make_shared<Eigen::VectorXd>();
};

// Obtain global nodal and pipe indices and pointers from meshes
void MatrixAssembler::global_nodal_pipe_indices(std::shared_ptr<Mesh> mesh) {
  global_nodes_.clear();
  global_pipes_.clear();
  for (const auto& node : mesh->nodes_) {
    global_nodes_.emplace(node);
  }
  for (const auto& pipe : mesh->pipes_) {
    global_pipes_.emplace(pipe);
  }
  nnode_ = global_nodes_.size();
  npipe_ = global_pipes_.size();
}

//! Assign pipe roughness coefficient for Hazen-Williams equation
void MatrixAssembler::assign_pipe_roughness(
    std::shared_ptr<Mesh> mesh,
    std::vector<std::pair<Index, double>> pipe_roughness) {
  for (auto& roughness : pipe_roughness) {
    mesh->pipes_.at(roughness.first)->pipe_roughness(roughness.second);
  }
}

//! Initialize discharges in pipes
void MatrixAssembler::initialize_pipe_discharge(std::shared_ptr<Mesh> mesh) {
  for (auto& pipe : mesh->pipes_) {
    pipe.second->initialize_discharge();
  }
}

//! Assign initial heads for nodes that have known head
void MatrixAssembler::assign_node_head(
    std::shared_ptr<Mesh> mesh,
    std::vector<std::pair<Index, double>> node_head) {
  for (auto& head : node_head) {
    mesh->nodes_.at(head.first)->head(head.second);
  }
}

//! Assign initial discharges for nodes that have known discharge
void MatrixAssembler::assign_node_discharge(
    std::shared_ptr<Mesh> mesh,
    std::vector<std::pair<Index, double>> node_discharge) {
  for (auto& discharge : node_discharge) {
    mesh->nodes_.at(discharge.first)->discharge(discharge.second);
  }
}

// Initialize nodal head vector
// If head of the ndoe is unknown (hasn't been assigned), initialize to zero
void MatrixAssembler::assemble_node_head_vector() {
  (*node_head_vec_).resize(nnode_);
  for (auto& node : global_nodes_) {
    Index index = node.first;
    if ((node.second)->ishead()) {
      (*node_head_vec_)(index) = node.second->head();
    } else {
      (*node_head_vec_)(index) = 0.0;
    }
  }
}

// Initialize nodal discharge vector
// If discharge of the ndoe is unknown (hasn't been assigned), initialize to
// zero
void MatrixAssembler::assemble_node_discharge_vector() {
  (*node_discharge_vec_).resize(nnode_);
  for (auto& node : global_nodes_) {
    Index index = node.first;
    if ((node.second)->isdischarge()) {
      (*node_discharge_vec_)(index) = node.second->discharge();
    } else {
      (*node_discharge_vec_)(index) = 0.0;
    }
  }
}

// Apply head to nodes
void MatrixAssembler::apply_node_head() {
  // Iterate through solved nodal head vector
  for (auto& node : global_nodes_) {
    // Get global index
    Index index = node.first;
    // Assign head
    node.second->head((*node_head_vec_)(index));
  }
}

// Apply discharge to nodes
void MatrixAssembler::apply_node_discharge() {
  // Iterate through solved nodal discharge vector
  for (auto& node : global_nodes_) {
    // Get global index
    Index index = node.first;
    // Assign discharge
    node.second->discharge((*node_discharge_vec_)(index));
  }
}

// Initialize pipe discharge vector
// Calculated according to nodal heads at two end and Hazen-Williams equation
// If any of the nodal head is unknown, initialize the discharge to 0.001
void MatrixAssembler::assemble_pipe_discharge_vector() {
  (*pipe_discharge_vec_).resize(npipe_);
  for (auto& pipe : global_pipes_) {
    Index index = pipe.first;
    (*pipe_discharge_vec_)(index) = 0.001;
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
//            start node, 0 for end node. 
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
void MatrixAssembler::assemble_jacobian() {
  // Check network
  if (nnode_ <= 0 || npipe_ <= 0) {
    throw std::runtime_error(
        "No node or pipe index pairs created to assemble Jacobian matrix");
  }

  std::vector<Eigen::Triplet<double>> update;
  update.reserve(nnode_ + 2 * npipe_);

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
    if (pipe.second->discharge() >= 0) {
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
