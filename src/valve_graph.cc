
#include "valve_graph.h"

pipenetwork::ValveGraph::ValveGraph(
    const std::shared_ptr<pipenetwork::MeshNodes>& mesh_nodes,
    const std::shared_ptr<pipenetwork::MeshLinks>& mesh_links,
    const std::vector<pipenetwork::isolation::ISOVProp>& iso_valve_props)
    : MeshGraph(mesh_nodes, mesh_links), iso_valves_(iso_valve_props) {
  construct_node_pipe_mtx();
  construct_valve_loc_mtx();
  construct_valve_def_mtx();
  initialize_iso_segs();
}

void pipenetwork::ValveGraph::construct_node_pipe_mtx() {
  std::vector<Eigen::Triplet<int>> graph_triplet;
  int connectivity_val = 1;

  for (const auto& index_link : mesh_links_->links()) {
    auto lid = index_link.first;
    auto link = index_link.second;
    auto end_nodes = link->nodes();
    graph_triplet.emplace_back(end_nodes.first->id(), lid, connectivity_val);
    graph_triplet.emplace_back(end_nodes.second->id(), lid, connectivity_val);
  }
  auto nnodes = mesh_nodes_->nnodes();
  auto nlinks = mesh_links_->nlinks();
  node_pipe_mtx_.resize(nnodes, nlinks);
  node_pipe_mtx_.setFromTriplets(graph_triplet.begin(), graph_triplet.end());
}

void pipenetwork::ValveGraph::construct_valve_loc_mtx() {
  auto iso_valve_props = iso_valves_.iso_valves();
  std::vector<Eigen::Triplet<int>> graph_triplet;
  int on_val = 1;
  for (const auto& vid_vprop : iso_valve_props) {
    auto vid = vid_vprop.first;
    auto vprop = vid_vprop.second;
    auto on_node = mesh_nodes_->get_node(vprop.on_node);
    auto on_link = mesh_links_->get_link(vprop.on_pipe);
    graph_triplet.emplace_back(on_node->id(), on_link->id(), on_val);

    iso_valves_.loc2vid(
        std::make_pair<Index, Index>(on_node->id(), on_link->id()), vid);
    mtx_helper.valve_loc_row2col[on_node->id()].emplace_back(on_link->id());
    mtx_helper.valve_loc_col2row[on_link->id()].emplace_back(on_node->id());
  }
  auto nnodes = mesh_nodes_->nnodes();
  auto nlinks = mesh_links_->nlinks();
  valve_loc_mtx_.resize(nnodes, nlinks);
  valve_loc_mtx_.setFromTriplets(graph_triplet.begin(), graph_triplet.end());
}

void pipenetwork::ValveGraph::construct_valve_def_mtx() {
  valve_def_mtx_ = node_pipe_mtx_ - valve_loc_mtx_;
  construct_vdef_idx_table();
}

void pipenetwork::ValveGraph::construct_vdef_idx_table() {
  valve_def_mtx_.makeCompressed();
  for (int k = 0; k < valve_def_mtx_.outerSize(); ++k)
    for (Eigen::SparseMatrix<int>::InnerIterator it(valve_def_mtx_, k); it;
         ++it) {
      auto row = it.row();  // row index
      auto col = it.col();  // col index (here it is equal to k)
      auto val = it.value();
      if (val != 0) {
        mtx_helper.valve_def_row2col[row].emplace_back(col);
        mtx_helper.valve_def_col2row[col].emplace_back(row);
      }
    }
}

void pipenetwork::ValveGraph::initialize_iso_segs() {
  std::set<Index> unexplored_pids;
  for (const auto& index_link : mesh_links_->links()) {
    auto lid = index_link.first;
    unexplored_pids.insert(lid);
  }
  iso_segments_ = isolation::IsoSegments();  // reinitialize for clear state
  iso_segments_.construct_iso_segs(mtx_helper, iso_valves_, unexplored_pids);
}

pipenetwork::isolation::IsoSeg pipenetwork::ValveGraph::pname2segment(
    const std::string& p_name) {
  auto pid = lname2id(p_name);
  return iso_segments_.pid2seg(pid);
}

std::vector<std::vector<pipenetwork::isolation::IsoSeg>>
    pipenetwork::ValveGraph::segment_components(
        const std::vector<std::string>& pipe2iso) {
  std::set<Index> sids;
  for (auto pname : pipe2iso) {
    auto seg = pname2segment(pname);
    sids.insert(seg.sid);
  }
  return iso_segments_.get_segment_components(sids);
}

void pipenetwork::ValveGraph::fail_valves(
    const std::vector<std::string>& broken_valves) {
  for (auto& vname : broken_valves) {
    auto vid = iso_valves_.vname2vid(vname);
    iso_segments_.merge_segments(vid);
  }
}
