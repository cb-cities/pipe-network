// Create nodal pointers and assign indices based on input coordinates
void Mesh::create_nodes(const std::vector<Eigen::Vector3d>& coords) {
  Index index = 0;
  for (const auto& coord : coords) {
    nodes_.emplace(std::pair<Index, std::shared_ptr<pipenetwork::Node>>(
        index, std::make_shared<pipenetwork::Node>(index, coord)));
    ++index;
  }
}

// Create pipe pointers and assign indices based on the nodes at its ends
bool Mesh::create_pipes(const std::vector<std::pair<Index, Index>>& nodeids,
                        const std::vector<double>& diameter,
                        const std::vector<double>& roughness,
                        const std::vector<bool>& pipe_status) {
  bool status = true;
  Index index = 0;
  try {
    for (const auto& nodeid : nodeids) {
      std::array<std::shared_ptr<pipenetwork::Node>, 2> nodes;
      nodes.at(0) = nodes_.at(nodeid.first);
      nodes.at(1) = nodes_.at(nodeid.second);
      pipes_.emplace(std::pair<Index, std::shared_ptr<pipenetwork::Pipe>>(
          index, std::make_shared<pipenetwork::Pipe>(
                     index, nodes, diameter.at(index), roughness.at(index),
                     pipe_status.at(index))));
      ++index;
    }
  } catch (std::out_of_range& range_error) {
    status = false;
    std::cout << "Pipe(id: " << index
              << ") is not created, as input node does not exist, nodal id: "
              << nodeids.at(index).first << " or " << nodeids.at(index).second
              << '\n'
              << "Pipe creation is unfinished" << '\n';
  }
  return status;
}

// Remove unconnected nodes from the mesh
void Mesh::remove_unconnected_nodes() {
  // Search for and record all unconnected nodes
  std::vector<Index> unconnected_nodes;
  for (const auto& node : nodes_) {
    bool is_connected = false;
    for (const auto& pipe : pipes_) {
      if (pipe.second->nodes().at(0)->id() == node.second->id() ||
          pipe.second->nodes().at(1)->id() == node.second->id()) {
        is_connected = true;
        break;
      }
    }
    if (!is_connected) unconnected_nodes.emplace_back(node.first);
  }
  // Remove isolated nodes from list of nodal indices
  for (const auto& unconnected_node : unconnected_nodes)
    nodes_.erase(unconnected_node);
}

// Initialize discharges in pipes
void Mesh::initialize_pipe_discharge() {
  for (auto& pipe : pipes_) pipe.second->initialize_discharge();
}

// Assign initial heads for nodes that have known head
void Mesh::assign_node_head(
    const std::vector<std::pair<Index, double>>& node_head) {
  for (auto& head : node_head) nodes_.at(head.first)->head(head.second);
}

// Assign initial discharges for nodes that have known discharge
void Mesh::assign_node_discharge(
    const std::vector<std::pair<Index, double>>& node_discharge) {
  for (auto& discharge : node_discharge)
    nodes_.at(discharge.first)->discharge(discharge.second);
}
