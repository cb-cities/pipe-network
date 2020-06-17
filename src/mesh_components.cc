#include "mesh_components.h"

pipenetwork::MeshNodes::MeshNodes(
    const std::vector<pipenetwork::JunctionProp>& junc_props,
    const std::vector<pipenetwork::ReservoirProp>& res_props) {

  add_nodes(junc_props);
  add_nodes(res_props);
}

std::shared_ptr<pipenetwork::Node> pipenetwork::MeshNodes::get_node(
    const std::string& node_name) const {
  try {
    Index nid = name2id_.at(node_name);
    return nodes_.at(nid);
  } catch (...) {
    throw std::runtime_error("Node does not exist: " + node_name + "\n");
  }
}
template <typename Prop>
void pipenetwork::MeshNodes::add_nodes(const std::vector<Prop>& props) {
  for (const auto& prop : props) {
    Index nid = nid_manager_.create_index();
    name2id_.emplace(prop.name, nid);
    add_node(prop);
  }
}

void pipenetwork::MeshNodes::add_node(
    const pipenetwork::JunctionProp& junc_prop) {
  auto nid = nid_manager_.current_index();
  auto junc = std::make_shared<pipenetwork::Junction>(nid, junc_prop);
  junctions_.emplace(nid, junc);
  nodes_.emplace(nid, junc);
}

void pipenetwork::MeshNodes::add_node(
    const pipenetwork::ReservoirProp& res_prop) {
  auto nid = nid_manager_.current_index();
  auto res = std::make_shared<pipenetwork::Reservoir>(nid, res_prop);
  reservoirs_.emplace(nid, res);
  nodes_.emplace(nid, res);
}

// Meshlinks constructor
pipenetwork::MeshLinks::MeshLinks(
    const std::vector<pipenetwork::PipeProp>& pipe_props,
    const std::vector<pipenetwork::PumpProp>& pump_props,
    const std::vector<pipenetwork::ValveProp>& valve_props,
    const MeshNodes& mesh_nodes) {

  add_links(pipe_props, mesh_nodes);
  add_links(pump_props, mesh_nodes);
  add_links(valve_props, mesh_nodes);
}

template <typename Prop>
void pipenetwork::MeshLinks::add_links(
    const std::vector<Prop>& props, const pipenetwork::MeshNodes& mesh_nodes) {
  for (const auto& prop : props) {
    try {
      auto node1 = mesh_nodes.get_node(prop.node1_name);
      auto node2 = mesh_nodes.get_node(prop.node2_name);
      if (prop.status != LinkStatus::CLOSED) {
        Index lid = lid_manager_.create_index();
        add_link(node1, node2, prop);
      }
    } catch (std::exception& e) {
      std::cout << "Failed to create link: " << prop.name << "\n";
      std::cout << "Exception: " << e.what() << "\n";
    }
  }
}

void pipenetwork::MeshLinks::add_link(
    const std::shared_ptr<pipenetwork::Node>& node1,
    const std::shared_ptr<pipenetwork::Node>& node2,
    const pipenetwork::PipeProp& pipe_prop) {
  auto lid = lid_manager_.current_index();
  auto pipe =
      std::make_shared<pipenetwork::Pipe>(lid, *node1, *node2, pipe_prop);
  pipes_.emplace(lid, pipe);
  links_.emplace(lid, pipe);
}

void pipenetwork::MeshLinks::add_link(
    const std::shared_ptr<pipenetwork::Node>& node1,
    const std::shared_ptr<pipenetwork::Node>& node2,
    const pipenetwork::PumpProp& pump_prop) {
  auto lid = lid_manager_.current_index();
  auto pump =
      std::make_shared<pipenetwork::Pump>(lid, *node1, *node2, pump_prop);
  pumps_.emplace(lid, pump);
  links_.emplace(lid, pump);
}

void pipenetwork::MeshLinks::add_link(
    const std::shared_ptr<pipenetwork::Node>& node1,
    const std::shared_ptr<pipenetwork::Node>& node2,
    const pipenetwork::ValveProp& valve_prop) {
  auto lid = lid_manager_.current_index();
  auto valve =
      std::make_shared<pipenetwork::Valve>(lid, *node1, *node2, valve_prop);
  valves_.emplace(lid, valve);
  links_.emplace(lid, valve);
}

void pipenetwork::MeshGraph::compute_graph_info_() {
  std::vector<Eigen::Triplet<int>> graph_triplet;
  int connectivity_val = 1;
  auto nnodes = mesh_nodes_->nnodes();
  for (const auto& index_link : mesh_links_->links()) {
    auto lid = index_link.first;
    auto link = index_link.second;
    auto end_nodes = link->nodes();
    node2link_[end_nodes.first->id()].emplace_back(lid);
    node2link_[end_nodes.second->id()].emplace_back(lid);

    graph_triplet.emplace_back(end_nodes.first->id(), end_nodes.second->id(),
                               connectivity_val);
    graph_triplet.emplace_back(end_nodes.second->id(), end_nodes.first->id(),
                               connectivity_val);
  }
  A_.resize(nnodes, nnodes);
  A_.setFromTriplets(graph_triplet.begin(), graph_triplet.end());

  compute_node_degrees_();
}

void pipenetwork::MeshGraph::compute_node_degrees_() {
  Index nnodes = mesh_nodes_->nnodes();
  ndegree_.resize(nnodes);
  for (int i = 0; i < nnodes; ++i) {
    ndegree_[i] = A_.outerIndexPtr()[i + 1] - A_.outerIndexPtr()[i];
  }
}

Eigen::VectorXd pipenetwork::MeshGraph::bfs(pipenetwork::Index nid) {
  Index nnodes = mesh_nodes_->nnodes();
  Eigen::VectorXd check_result(nnodes);
  check_result.setZero(nnodes);

  std::set<int> nodes_to_explore;
  nodes_to_explore.emplace(nid);
  check_result[nid] = 1;
  auto indptr = A_.outerIndexPtr();
  auto indices = A_.innerIndexPtr();
  auto data = A_.valuePtr();

  while (!nodes_to_explore.empty()) {
    auto node_being_explored = *nodes_to_explore.begin();
    nodes_to_explore.erase(nodes_to_explore.begin());

    int nconnections = ndegree_[node_being_explored];
    int ndx = indptr[node_being_explored];
    // for all the connected nodes, set result to 1 and place them into the
    // searching queue
    for (int i = 0; i < nconnections; ++i) {
      if (data[ndx + i] != 0 && check_result[indices[ndx + i]] == 0) {
        check_result[indices[ndx + i]] = 1;
        nodes_to_explore.emplace(indices[ndx + i]);
      }
    }
  }
  return check_result;
}