//! Create nodal pointers and assign indices based on input coordinates
void Mesh::create_nodes(std::vector<Eigen::Vector3d> coords) {
  unsigned index = 1;
  for (auto& coord : coords) {
    nodes_.emplace(std::pair<unsigned, std::shared_ptr<pipenetwork::Node>>(
        index, std::make_shared<pipenetwork::Node>(index, coord)));
    index++;
  }
}

//! Create pipe pointers and assign indices based on the nodes at its ends
void Mesh::create_pipes(std::vector<std::pair<unsigned, unsigned>> nodeids) {
  unsigned index = 1;
  for (auto& nodeid : nodeids) {
    std::array<std::shared_ptr<pipenetwork::Node>, 2> nodes;
    nodes.at(0) = nodes_.at(nodeid.first);
    nodes.at(1) = nodes_.at(nodeid.second);
    pipes_.emplace(std::pair<unsigned, std::unique_ptr<pipenetwork::Pipe>>(
        index, std::make_unique<pipenetwork::Pipe>(index, nodes)));
    index++;
  }
}

//! Check whether isolated node exists
//! \retval the status to indicate whether isolated node exists
bool Mesh::isolated_node() {
  for (const auto& node : nodes_) {
    bool connect = false;
    for (const auto& pipe : pipes_) {
      if (pipe.second->nodes().at(0)->id() == node.second->id() ||
          pipe.second->nodes().at(1)->id() == node.second->id()) {
        connect = true;
        break;
      }
    }
    if (connect == false) return true;
  }
  return false;
}

//! Return coordinates of all the nodes in the mesh
//! \retval nodal_coordinates coordinates of all the nodes
std::vector<Eigen::Vector3d> Mesh::nodal_coordinates() {
  std::vector<Eigen::Vector3d> nodal_coordinates;
  for (const auto& node : nodes_)
    nodal_coordinates.emplace_back(node.second->coordinates());
  return nodal_coordinates;
}
