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

void pipenetwork::Mesh::find_iso_components_() {
  find_iso_nodes_();
  find_iso_links_();
}

void pipenetwork::Mesh::find_iso_nodes_() {
  Index nnodes = mesh_nodes_->nnodes();
  Eigen::VectorXd check_result(nnodes);
  check_result.setZero(nnodes);

  auto reservoirs = mesh_nodes_->reservoirs();
  for (const auto& index_res : reservoirs) {
    auto res_id = index_res.first;
    check_result += mesh_graph_->bfs(res_id);
  }
  for (Index i = 0; i < nnodes; ++i) {
    if (check_result[i] == 0) {
      iso_nodes_.emplace_back(i);
    }
  }
}

void pipenetwork::Mesh::find_iso_links_() {
  auto node2link_map = mesh_graph_->node2link_map();
  for (const auto& nid : iso_nodes_) {
    if (node2link_map.find(nid) != node2link_map.end()) {
      auto links = node2link_map.at(nid);
      for (const auto& link : links) {
        iso_links_.emplace_back(link);
      }
    }
  }
}
void pipenetwork::Mesh::print_summary() {
  std::cout << "Network Name: " << name_ << std::endl
            << "number of pipes: " << mesh_links_->npipes()
            << " ;number of pumps: " << mesh_links_->npumps()
            << " ;number of valves: " << mesh_links_->nvalves()
            << " ;number of junctions: " << mesh_nodes_->njunctions()
            << " ;number of sources: " << mesh_nodes_->nreservoirs()
            << std::endl;

  //  std::cout
  //      << "======================= NODES INFO ======================="
  //      << std::endl;
  //
  //  for (auto const& x : nodes_) {
  //    std::cout << "node id: " << x.second->id() << std::endl;
  //  }
  //
  //  std::cout << "======================= LINKS INFO
  //  == == == == == == == == == == == = "
  //            << std::endl;
  //  for (auto const& x : links_) {
  //    std::cout << "link id: " << x->id()
  //              << "; node1 id: " << x->nodes().first->id()
  //              << "; node2 id: " << x->nodes().second->id() << std::endl;
  //  }
}
