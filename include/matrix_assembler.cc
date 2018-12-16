// Constructor
MatrixAssembler::MatrixAssembler() {
  jac_ = std::make_shared<Eigen::SparseMatrix<double>>();
  node_head_vec_ = std::make_shared<Eigen::VectorXd>();
  node_discharge_vec_ = std::make_shared<Eigen::VectorXd>();
  pipe_discharge_vec_ = std::make_shared<Eigen::VectorXd>();
};

// Assign global index to nodes in the whole network
void MatrixAssembler::assign_global_nodal_indices(
    std::vector<std::shared_ptr<pipenetwork::Node>> nodes) {
  global_nodes_.clear();
  Index index = 0;
  for (const auto& node : nodes) {
    global_nodes_.emplace(
        std::pair<std::shared_ptr<pipenetwork::Node>, Index>(node, index));
    ++index;
  }
  nnode_ = global_nodes_.size();
}

// Assign global index to pipes in the whole network
void MatrixAssembler::assign_global_pipe_indices(
    std::vector<std::shared_ptr<pipenetwork::Pipe>> pipes) {
  global_pipes_.clear();
  Index index = 0;
  for (const auto& pipe : pipes) {
    global_pipes_.emplace(
        std::pair<std::shared_ptr<pipenetwork::Pipe>, Index>(pipe, index));
    ++index;
  }
  npipe_ = global_pipes_.size();
}

// Initialize nodal head vector
// If head of the ndoe is unknown (hasn't been assigned), initialize to zero
void MatrixAssembler::assemble_node_head_vector() {
  (*node_head_vec_).resize(nnode_);
  for (auto node : global_nodes_) {
    Index index = node.second;
    if ((node.first)->ishead()) {
      (*node_head_vec_)(index) = node.first->head();
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
  for (auto node : global_nodes_) {
    Index index = node.second;
    if ((node.first)->isdischarge()) {
      (*node_discharge_vec_)(index) = node.first->discharge();
    } else {
      (*node_discharge_vec_)(index) = 0.0;
    }
  }
}

// Apply head to nodes
void MatrixAssembler::apply_node_head() {
  // Iterate through solved nodal head vector
  for (auto node : global_nodes_) {
    // Get global index
    Index index = node.second;
    // Assign head
    node.first->head((*node_head_vec_)(index));
  }
}

// Apply discharge to nodes
void MatrixAssembler::apply_node_discharge() {
  // Iterate through solved nodal discharge vector
  for (auto node : global_nodes_) {
    // Get global index
    Index index = node.second;
    // Assign discharge
    node.first->discharge((*node_discharge_vec_)(index));
  }
}

// Initialize pipe discharge vector
// Calculated according to nodal heads at two end and Hazen-Williams equation
// If any of the nodal head is unknown, initialize the discharge to 0.001
void MatrixAssembler::assemble_pipe_discharge_vector() {
  (*pipe_discharge_vec_).resize(npipe_);
  for (auto pipe : global_pipes_) {
    Index index = pipe.second;
    (*pipe_discharge_vec_)(index) = 0.001;
  }
}

// Assemble Jacobian matrix
//                 nodal_head    nodal_discharge    pipe_discharge
// nadal_balance   sub_jac_A        sub_jac_B         sub_jac_C
// headloss        sub_jac_D        sub_jac_E         sub_jac_F
void MatrixAssembler::assemble_jacobian() {
  // Check network
  if (nnode_ <= 0 || npipe_ <= 0) {
    std::cerr
        << "Error: create node and pipe index pair to assemble Jacobian matrix"
        << std::endl;
  }

  std::vector<Eigen::Triplet<double>> update;
  update.reserve(nnode_ + 2 * npipe_);

  // Iterate through all nodes
  for (const auto& node : global_nodes_) {
    Index index = node.second;
    // construct jacB part
    update.emplace_back(index, nnode_ + index, -1);
  }

  // Iterate through all pipes
  for (const auto& pipe : global_pipes_) {
    Index index_pipe = pipe.second;
    Index index_node1 = 0;
    Index index_node2 = 0;
    for (const auto& node : global_nodes_) {
      if (node.first == pipe.first->nodes().at(0)) {
        index_node1 = node.second;
      }
      if (node.first == pipe.first->nodes().at(1)) {
        index_node2 = node.second;
      }
    }
    // Index index_node1 = global_nodes_.at(pipe.first->nodes().at(0));
    // Index index_node2 = global_nodes_.at(pipe.first->nodes().at(1));
    // construct jacC part
    if (pipe.first->discharge() >= 0) {
      update.emplace_back(index_node1, 2 * nnode_ + index_pipe, -1);
      update.emplace_back(index_node2, 2 * nnode_ + index_pipe, 1);
    } else {
      update.emplace_back(index_node1, 2 * nnode_ + index_pipe, -1);
      update.emplace_back(index_node2, 2 * nnode_ + index_pipe, 1);
    }
    // construct jacD part
    update.emplace_back(nnode_ + index_pipe, index_node1, 1);
    update.emplace_back(nnode_ + index_pipe, index_node2, -1);
    // construct jacF part (Hazen-Williams)
    double discharge = pipe.first->discharge();
    double coeff = 10.67 * pipe.first->length() /
                   (pow(pipe.first->pipe_roughness(), 1.852) *
                    pow(pipe.first->diameter(), 4.8704));
    double residual_deriv = -1.852 * coeff * pow(discharge, 0.852);
    update.emplace_back(nnode_ + index_pipe, 2 * nnode_ + index_pipe,
                        -1.852 * coeff * pow(discharge, 0.852));
  }

  jac_->resize(nnode_ + npipe_, 2 * nnode_ + npipe_);
  jac_->setFromTriplets(update.begin(), update.end());
}
