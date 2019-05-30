
#include "mesh.h"

// Create nodal pointers and assign indices based on input coordinates
void pipenetwork::Mesh::create_nodes(
    const std::vector<Eigen::Vector3d>& coords) {
  Index index = 0;
  for (const auto& coord : coords) {
    nodes_.emplace(std::pair<Index, std::shared_ptr<pipenetwork::Node>>(
        index, std::make_shared<pipenetwork::Node>(index, coord)));
    ++index;
  }
}

// Create nodal pointers and assign indices based on input ID and coordinates
void pipenetwork::Mesh::create_nodes(
    const std::vector<std::pair<Index, Eigen::Vector3d>>& coords_info) {
  for (const auto& coord_info : coords_info) {
    auto index = coord_info.first;
    auto coord = coord_info.second;
    nodes_.emplace(std::pair<Index, std::shared_ptr<pipenetwork::Node>>(
        index, std::make_shared<pipenetwork::Node>(index, coord)));
  }
}

// Create pipe pointers and assign indices based on the nodes at its ends
bool pipenetwork::Mesh::create_pipes(
    const std::vector<std::pair<Index, Index>>& nodeids,
    const std::vector<double>& diameter, const std::vector<double>& length,
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
                     index, nodes, diameter.at(index),length.at(index), roughness.at(index),
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

// Create pipe pointers and assign indices based on the nodes at its ends
bool pipenetwork::Mesh::create_pipes(
    const std::vector<std::tuple<int, int, int>>& ids,
    const std::vector<double>& diameter, const std::vector<double>& length,
    const std::vector<double>& roughness,
    const std::vector<bool>& pipe_status) {
  bool status = true;
  Index index = 0;
  try {
    for (const auto& id : ids) {

      std::array<std::shared_ptr<pipenetwork::Node>, 2> nodes;
      auto pipe_id = std::get<0>(id);
      nodes.at(0) = nodes_.at(std::get<1>(id));
      nodes.at(1) = nodes_.at(std::get<2>(id));
      pipes_.emplace(std::pair<Index, std::shared_ptr<pipenetwork::Pipe>>(
          pipe_id, std::make_shared<pipenetwork::Pipe>(
                       index, nodes, diameter.at(index), length.at(index), roughness.at(index),
                       pipe_status.at(index))));
      ++index;
    }
  } catch (std::out_of_range& range_error) {
    status = false;
    std::cout << "Pipe(vec loc: " << index
              << ") is not created, as input node does not exist" << '\n'
              << "Pipe creation is unfinished" << '\n';
  }
  return status;
}

// Remove unconnected nodes from the mesh
void pipenetwork::Mesh::remove_unconnected_nodes() {
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
void pipenetwork::Mesh::initialize_pipe_discharge(
    const std::vector<std::pair<Index, double>>& init_discharge) {
  for (auto& pipe : pipes_) pipe.second->initialize_discharge();
  for (const auto& discharge : init_discharge)
    pipes_.at(discharge.first)->initialize_discharge(discharge.second);
}

void pipenetwork::Mesh::initialize_pipe_discharge(double init_discharge) {
  for (auto& pipe : pipes_) pipe.second->initialize_discharge(init_discharge);
}

// Assign initial heads for nodes that have known head
void pipenetwork::Mesh::assign_node_head(
    const std::vector<std::pair<Index, double>>& node_head) {
  for (const auto& head : node_head) nodes_.at(head.first)->head(head.second);
}

// Assign initial elevation for nodes that have known head
void pipenetwork::Mesh::assign_node_elevation(
    const std::vector<std::pair<Index, double>>& node_head) {
  for (const auto& head : node_head)
    nodes_.at(head.first)->elevation(head.second);
}

// Assign initial discharges for nodes that have known discharge
void pipenetwork::Mesh::assign_node_demand(
    const std::vector<std::pair<Index, double>>& node_discharge) {
  for (const auto& discharge : node_discharge)
    nodes_.at(discharge.first)->demand(discharge.second);
}

// Assign leak diameter for nodes that have known leaks
void pipenetwork::Mesh::assign_node_leak(
    const std::vector<std::pair<Index, double>>& node_leak_dia) {
  for (const auto& leak : node_leak_dia)
    nodes_.at(leak.first)->leak_diameter(leak.second);
}