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
            << " ;number of pumps: " << npumps()
            << " ;number of valves: "<< nvalves()<< std::endl;

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

void pipenetwork::Mesh::create_valve(
    std::vector<pipenetwork::Valve_prop>& valve_props) {
  for (auto& valve_prop : valve_props) {
    auto node1id = valve_prop.node1_id;
    auto node2id = valve_prop.node2_id;

    connected_nodes_.emplace(node1id, nodes_.at(node1id));
    connected_nodes_.emplace(node2id, nodes_.at(node2id));

    valve_prop.node1 = nodes_.at(node1id);
    valve_prop.node2 = nodes_.at(node2id);

    links_.emplace_back(std::make_shared<pipenetwork::Valve>(valve_prop));
    ++nvalves_;
  }
}

void pipenetwork::Mesh::create_mesh_from_inp(
    std::shared_ptr<pipenetwork::Input>& IO) {
  auto pipe_props = IO->pipe_properties();
  auto pump_props = IO->pump_properties();
  auto valve_props = IO->valve_properties();

  // create nodes
  create_junctions(IO->junction_properties());
  create_reservoirs(IO->reservoir_properties());

  // create links
  create_pipes(pipe_props);
  if (!pump_props.empty()) {
    create_pumps(pump_props);
  }
  if (!valve_props.empty()) {
    create_valve (valve_props);
  }
}
