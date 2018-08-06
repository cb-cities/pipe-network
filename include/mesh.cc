//! Create nodal pointers and assign indices based on input coordinates
void Mesh::create_nodes(const std::vector<Eigen::Vector3d>& coords) {
  Index index = 0;
  for (const auto& coord : coords) {
    nodes_.emplace(std::pair<Index, std::shared_ptr<pipenetwork::Node>>(
        index, std::make_shared<pipenetwork::Node>(index, coord)));
    ++index;
  }
}

//! Create pipe pointers and assign indices based on the nodes at its ends
void Mesh::create_pipes(const std::vector<std::pair<Index, Index>>& nodeids) {
  Index index = 0;
  for (const auto& nodeid : nodeids) {
    std::array<std::shared_ptr<pipenetwork::Node>, 2> nodes;
    nodes.at(0) = nodes_.at(nodeid.first);
    nodes.at(1) = nodes_.at(nodeid.second);
    pipes_.emplace(std::pair<Index, std::unique_ptr<pipenetwork::Pipe>>(
        index, std::make_unique<pipenetwork::Pipe>(index, nodes)));
    ++index;
  }
}

//! Record isolated nodes and remove them from the mesh
std::vector<std::shared_ptr<pipenetwork::Node>> Mesh::isolated_nodes() {
  // Search for and collect all isolated nodes
  std::vector<std::shared_ptr<pipenetwork::Node>> isolated_nodes;
  for (const auto& node : nodes_) {
    bool connect = false;
    for (const auto& pipe : pipes_) {
      if (pipe.second->nodes().at(0)->id() == node.second->id() ||
          pipe.second->nodes().at(1)->id() == node.second->id()) {
        connect = true;
        break;
      }
    }
    if (connect == false) isolated_nodes.emplace_back(node.second);
  }
  // Remove isolated nodes from list of nodal pointers
  for (const auto& inode : isolated_nodes) nodes_.erase(inode->id());
  // Return the list of isolated nodes
  return isolated_nodes;
}
