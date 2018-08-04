//! Create a nodal pointer based on input id and coordinates and check whether
//! duplicate nodes exist
void Mesh::create_node(std::vector<unsigned> id,
                       std::vector<Eigen::Vector3d> coords) {
  // Check whether there are duplicate indices or coordinates
  if (id.size() > 1) {
    for (auto iter1 = id.begin(); iter1 != std::prev(id.end()); ++iter1) {
      for (auto iter2 = std::next(iter1); iter2 != id.end(); ++iter2) {
        if (*iter1 == *iter2) {
          std::cout << "node id: " << *iter1 << '\n';
          throw std::runtime_error(
              "duplicate nodal indices found, check input");
        }
      }
    }
  }
  if (coords.size() > 1) {
    unsigned count = 0;
    for (auto iter1 = coords.begin(); iter1 != std::prev(coords.end());
         ++iter1) {
      for (auto iter2 = std::next(iter1); iter2 != coords.end(); ++iter2) {
        if (*iter1 == *iter2) {
          std::cout << "node id: " << id.at(count) << '\n';
          throw std::runtime_error(
              "duplicate nodal coordinates found, check input");
        }
      }
      count++;
    }
  }
  // create nodal pointer and put nodal id and corresponding pointer into map
  for (unsigned i = 0; i < id.size(); ++i)
    nodes_.insert(std::pair<unsigned, std::shared_ptr<pipenetwork::Node>>(
        id.at(i), std::make_shared<pipenetwork::Node>(id.at(i), coords.at(i))));
}

//! Create a pipe pointer based on input idices of the pipe and the nodes at its
//! ends and check whether duplicate pipe exist
void Mesh::create_pipe(std::vector<unsigned> pipeid,
                       std::vector<std::pair<unsigned, unsigned>> nodeids) {
  // Check whether there are duplicate indices
  if (pipeid.size() > 1) {
    for (auto iter1 = pipeid.begin(); iter1 != std::prev(pipeid.end());
         ++iter1) {
      for (auto iter2 = std::next(iter1); iter2 != pipeid.end(); ++iter2) {
        if (*iter1 == *iter2) {
          std::cout << "pipe id: " << *iter1 << '\n';
          throw std::runtime_error("duplicate pipe indices found, check input");
        }
      }
    }
  }
  // create pipe pointer and put pipe id and corresponding pointer into map
  for (unsigned i = 0; i < pipeid.size(); ++i) {
    std::array<std::shared_ptr<pipenetwork::Node>, 2> nodes;
    nodes.at(0) = nodes_.at(nodeids.at(i).first);
    nodes.at(1) = nodes_.at(nodeids.at(i).second);
    pipes_.insert(std::pair<unsigned, std::unique_ptr<pipenetwork::Pipe>>(
        pipeid.at(i),
        std::make_unique<pipenetwork::Pipe>(pipeid.at(i), nodes)));
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
