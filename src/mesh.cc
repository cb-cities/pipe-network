#include "mesh.h"

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
