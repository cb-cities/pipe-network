#include "mesh.h"

//
// void pipenetwork::Mesh::create_pipes(std::vector<Pipe_prop>& pipe_props) {
//  for (auto& pipe_prop : pipe_props) {
//    auto node1id = pipe_prop.node1_id;
//    auto node2id = pipe_prop.node2_id;
//
//    pipe_prop.node1 = nodes_.at(node1id);
//    pipe_prop.node2 = nodes_.at(node2id);
//
//    links_.emplace_back(std::make_shared<pipenetwork::Pipe>(pipe_prop));
//    ++npipes_;
//  }
//}
//
// void pipenetwork::Mesh::print_summary() {
//
//  std::cout << "number of nodes: " << nnodes()
//            << " ;number of links: " << nlinks()
//            << " ;number of pipes: " << npipes()
//            << " ;number of pumps: " << npumps()
//            << " ;number of valves: " << nvalves()
//            << " ;number of junctions: " << njunctions_
//            << " ;number of sources: " << nsrcs_ << std::endl;
//
//  //  std::cout
//  //      << "======================= NODES INFO ======================="
//  //      << std::endl;
//  //
//  //  for (auto const& x : nodes_) {
//  //    std::cout << "node id: " << x.second->id() << std::endl;
//  //  }
//  //
//  //  std::cout << "======================= LINKS INFO
//  ======================="
//  //            << std::endl;
//  //  for (auto const& x : links_) {
//  //    std::cout << "link id: " << x->id()
//  //              << "; node1 id: " << x->nodes().first->id()
//  //              << "; node2 id: " << x->nodes().second->id() << std::endl;
//  //  }
//}
//
// void pipenetwork::Mesh::create_pumps(
//    std::vector<pipenetwork::Pump_prop>& pump_props) {
//  for (auto& pump_prop : pump_props) {
//    auto node1id = pump_prop.node1_id;
//    auto node2id = pump_prop.node2_id;
//
//    pump_prop.node1 = nodes_.at(node1id);
//    pump_prop.node2 = nodes_.at(node2id);
//
//    links_.emplace_back(std::make_shared<pipenetwork::Pump>(pump_prop));
//    ++npumps_;
//  }
//}
//
// void pipenetwork::Mesh::create_valve(
//    std::vector<pipenetwork::Valve_prop>& valve_props) {
//  for (auto& valve_prop : valve_props) {
//    auto node1id = valve_prop.node1_id;
//    auto node2id = valve_prop.node2_id;
//
//    valve_prop.node1 = nodes_.at(node1id);
//    valve_prop.node2 = nodes_.at(node2id);
//
//    links_.emplace_back(std::make_shared<pipenetwork::Valve>(valve_prop));
//    ++nvalves_;
//  }
//}
//
// void pipenetwork::Mesh::create_mesh_from_inp(
//    std::shared_ptr<pipenetwork::Input>& IO) {
//  auto pipe_props = IO->pipe_properties();
//  auto pump_props = IO->pump_properties();
//  auto valve_props = IO->valve_properties();
//
//  // create nodes
//  try {
//    create_junctions(IO->junction_properties());
//  } catch (std::exception& e) {
//    std::cerr
//        << "Failed to create junctions from the input file, error message "
//        << e.what() << std::endl;
//    std::abort();
//  }
//  try {
//    create_reservoirs(IO->reservoir_properties());
//  } catch (std::exception& e) {
//    std::cerr
//        << "Failed to create reservoirs from the input file, error message "
//        << e.what() << std::endl;
//    std::abort();
//  }
//  // create links
//  try {
//    create_pipes(pipe_props);
//  } catch (std::exception& e) {
//    std::cerr << "Failed to create pipes from the input file, error message "
//              << e.what() << std::endl;
//    std::abort();
//  }
//
//  if (!pump_props.empty()) {
//    try {
//      create_pumps(pump_props);
//    } catch (std::exception& e) {
//      std::cerr << "Failed to create pumps from the input file, error message
//      "
//                << e.what() << std::endl;
//      std::abort();
//    }
//  }
//  if (!valve_props.empty()) {
//    try {
//      create_valve(valve_props);
//    } catch (std::exception& e) {
//      std::cerr << "Failed to create valves from the input file, error message
//      "
//                << e.what() << std::endl;
//      std::abort();
//    }
//  }
//}

pipenetwork::MeshNodes::MeshNodes(
    const std::vector<pipenetwork::JunctionProp>& junc_props,
    const std::vector<pipenetwork::ReservoirProp>& res_props) {

  add_nodes(junc_props);
  add_nodes(res_props);
}

std::shared_ptr<pipenetwork::Node> pipenetwork::MeshNodes::get_node(
    const std::string& node_name) const {
  try {
    Index nid = name2id.at(node_name);
    if (nid < junctions.size()) {
      return junctions.at(nid);
    }
    return reservoirs.at(nid);

  } catch (...) {
    throw std::runtime_error("Node does not exist: " + node_name + "\n");
  }
}
template <typename Prop>
void pipenetwork::MeshNodes::add_nodes(const std::vector<Prop>& props) {
  for (const auto& prop : props) {
    Index nid = nid_manager_.create_index();
    name2id.emplace(prop.name, nid);
    add_node(prop);
  }
}

void pipenetwork::MeshNodes::add_node(
    const pipenetwork::JunctionProp& junc_prop) {
  auto nid = nid_manager_.current_index();
  auto junc = std::make_shared<pipenetwork::Junction>(nid, junc_prop);
  junctions.emplace(nid, junc);
}

void pipenetwork::MeshNodes::add_node(
    const pipenetwork::ReservoirProp& res_prop) {
  auto nid = nid_manager_.current_index();
  auto res = std::make_shared<pipenetwork::Reservoir>(nid, res_prop);
  reservoirs.emplace(nid, res);
}

// Meshlinks constructor
pipenetwork::MeshLinks::MeshLinks(
    std::vector<pipenetwork::PipeProp>& pipe_props,
    std::vector<pipenetwork::PumpProp>& pump_props,
    std::vector<pipenetwork::ValveProp>& valve_props,
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
      Index lid = lid_manager_.create_index();
      auto node1 = mesh_nodes.get_node(prop.node1_name);
      auto node2 = mesh_nodes.get_node(prop.node2_name);
      add_link(node1, node2, prop);
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
  pipes.emplace(lid, pipe);
}

void pipenetwork::MeshLinks::add_link(
    const std::shared_ptr<pipenetwork::Node>& node1,
    const std::shared_ptr<pipenetwork::Node>& node2,
    const pipenetwork::PumpProp& pump_prop) {
  auto lid = lid_manager_.current_index();
  auto pump =
      std::make_shared<pipenetwork::Pump>(lid, *node1, *node2, pump_prop);
  pumps.emplace(lid, pump);
}

void pipenetwork::MeshLinks::add_link(
    const std::shared_ptr<pipenetwork::Node>& node1,
    const std::shared_ptr<pipenetwork::Node>& node2,
    const pipenetwork::ValveProp& valve_prop) {
  auto lid = lid_manager_.current_index();
  auto valve =
      std::make_shared<pipenetwork::Valve>(lid, *node1, *node2, valve_prop);
  valves.emplace(lid, valve);
}


