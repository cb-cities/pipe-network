#include "mesh.h"

void pipenetwork::Mesh::create_junctions(
    const std::vector<Junction_prop>& junc_props) {
  for (const auto& junc_prop : junc_props) {
    nodes_.emplace(junc_prop.id,
                   std::make_shared<pipenetwork::Junction>(junc_prop));
  }
}

void pipenetwork::Mesh::create_reservoirs(
    const std::vector<Reservoir_prop>& res_props) {
  for (const auto& res_prop : res_props) {
    nodes_.emplace(res_prop.id,
                   std::make_shared<pipenetwork::Reservoir>(res_prop));
  }
}

void pipenetwork::Mesh::create_pipes(std::vector<Pipe_prop>& pipe_props) {
  for (auto& pipe_prop : pipe_props) {
    auto node1id = pipe_prop.node1_id;
    auto node2id = pipe_prop.node2_id;
    connected_nodes_.emplace(node1id, nodes_.at(node1id));
    connected_nodes_.emplace(node2id, nodes_.at(node2id));

    pipe_prop.node1 = nodes_.at(node1id);
    pipe_prop.node2 = nodes_.at(node2id);

    links_.emplace_back(std::make_shared<pipenetwork::Pipe>(pipe_prop));
    ++npipes_;
  }
}

void pipenetwork::Mesh::print_summary() {
  auto nnode = nodes_.size();

  std::cout << "number of nodes: " << nnode << " ;number of links: " << nlinks()
            << " ;number of connected nodes: " << nnodes()
            << " ;number of pipes: " << npipes()
            << " ;number of pumps: " << npumps() << std::endl;

  std::cout
      << "======================= CONNECTED NODE INFO ======================="
      << std::endl;

  for (auto const& x : connected_nodes_) {
    std::cout << "node id: " << x.second->id() << std::endl;
  }

  std::cout << "======================= LINK INFO ======================="
            << std::endl;
  for (auto const& x : links_) {
    std::cout << "link id: " << x->id()
              << "; node1 id: " << x->nodes().first->id()
              << "; node2 id: " << x->nodes().second->id() << std::endl;
  }
}

void pipenetwork::Mesh::create_pumps(
    std::vector<pipenetwork::Pump_prop>& pump_props) {
  for (auto& pump_prop : pump_props) {
    auto node1id = pump_prop.node1_id;
    auto node2id = pump_prop.node2_id;
    connected_nodes_.emplace(node1id, nodes_.at(node1id));
    connected_nodes_.emplace(node2id, nodes_.at(node2id));

    pump_prop.node1 = nodes_.at(node1id);
    pump_prop.node2 = nodes_.at(node2id);

    links_.emplace_back(std::make_shared<pipenetwork::Pump>(pump_prop));
    ++npumps_;
  }
}
