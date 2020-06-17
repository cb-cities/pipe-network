#include "mesh.h"

void pipenetwork::Mesh::create_mesh_from_io(
    const std::shared_ptr<pipenetwork::IO>& IO) {
  create_nodes(IO->junction_properties(), IO->reservoir_properties());
  create_links(IO->pipe_properties(), IO->pump_properties(),
               IO->valve_properties());
  create_mesh_graph();
}

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

void pipenetwork::Mesh::find_leak_nids() {
  auto junctions = mesh_nodes_->junctions();
  for (const auto& index_junc : junctions) {
    auto nid = index_junc.first;
    auto junction = index_junc.second;
    if (junction->leak_area() > 0) {
      leak_nids_.emplace_back(nid);
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
            << " ;number of leaking junctions: " << leak_nids_.size()
            << " ;number of isolated junctions: " << iso_nodes_.size()
            << " ;number of isolated links: " << iso_links_.size() << std::endl;
}

void pipenetwork::Mesh::save_mesh(const std::string& output_path) {
  std::ofstream outnode(output_path + name_ + "_nodes.csv");
  std::ofstream outlink(output_path + name_ + "_links.csv");
  outnode << "node_name"
          << ","
          << "head"
          << ","
          << "demand"
          << "\n";
  outlink << "link_name"
          << ","
          << "flowrate"
          << "\n";

  // junctions
  auto junction_map = mesh_nodes_->junctions();
  for (auto& index_junc : junction_map) {
    auto nid = index_junc.first;
    auto junction = index_junc.second;
    auto junc_prop = junction->property();
    outnode << std::setprecision(12) << junc_prop.name << "," << junction->head
            << "," << junction->demand << "\n";
  }

  auto pipe_map = mesh_links_->pipes();
  for (auto& index_pipe : pipe_map) {
    auto nid = index_pipe.first;
    auto pipe = index_pipe.second;
    auto pipe_prop = pipe->property();
    outlink << std::setprecision(12) << pipe_prop.name << "," << pipe->flowrate
            << "\n";
  }
}
