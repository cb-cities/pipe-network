
#include "valve_graph.h"

pipenetwork::ValveGraph::ValveGraph(
    const std::shared_ptr<pipenetwork::MeshNodes>& mesh_nodes,
    const std::shared_ptr<pipenetwork::MeshLinks>& mesh_links,
    const std::vector<pipenetwork::ISOVProp>& iso_valve_props)
    : MeshGraph(mesh_nodes, mesh_links) {
  construct_node_pipe_mtx();
  construct_valve_loc_mtx(iso_valve_props);
  construct_valve_def_mtx();
  construct_iso_segs();
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

void pipenetwork::ValveGraph::construct_valve_loc_mtx(
    const std::vector<pipenetwork::ISOVProp>& iso_valve_props) {
  std::vector<Eigen::Triplet<int>> graph_triplet;
  int on_val = 1;
  for (const auto& vprop : iso_valve_props) {
    auto on_node = mesh_nodes_->get_node(vprop.on_node);
    auto on_link = mesh_links_->get_link(vprop.on_pipe);
    graph_triplet.emplace_back(on_node->id(), on_link->id(), on_val);
    mtx_helper
        .id2vname[std::make_pair<Index, Index>(on_node->id(), on_link->id())] =
        vprop.name;
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

pipenetwork::IsoSeg pipenetwork::ValveGraph::get_single_seg(
    pipenetwork::Index pid) {
  IsoSeg seg;
  bool iso = mtx_helper.valve_def_col2row.find(pid) ==
             mtx_helper.valve_def_col2row.end();
  if (iso) {
    seg.pids.insert(pid);
    return seg;
  }
  std::set<Index> pids{pid};
  while (!pids.empty()) {
    // col search (find isolated nodes)
    auto searching_pid = *pids.begin();
    pids.erase(pids.begin());
    seg.pids.insert(searching_pid);
    // adding new nids
    auto new_nids = mtx_helper.valve_def_col2row[pid];
    // row search (find new isolated pipes)
    while (!new_nids.empty()) {
      auto searching_nid = new_nids.back();
      new_nids.pop_back();
      seg.nids.insert(searching_nid);
      // adding new pids
      auto new_pids = mtx_helper.valve_def_row2col[searching_nid];
      for (auto new_pid : new_pids) {
        bool psearched = std::find(seg.pids.begin(), seg.pids.end(), new_pid) !=
                         seg.pids.end();
        if (!psearched) {
          pids.insert(new_pid);
        }
      }
    }
  }
  seg.vnames = get_seg_valves(seg.pids, seg.nids);
  return seg;
}

std::set<std::string> pipenetwork::ValveGraph::get_seg_valves(
    std::set<Index>& pids, std::set<Index>& nids) {
  std::set<std::string> valves;
  // pick by pipes
  for (Index pid : pids) {
    auto nids = mtx_helper.valve_loc_col2row[pid];
    for (Index nid : nids) {
      auto vid = std::make_pair<long, long>(nid, pid);
      auto valve = mtx_helper.id2vname[vid];
      //      std::cout<<"p pick"<<valve<<std::endl;
      valves.insert(valve);
    }
  }
  // pick by nodes
  for (Index nid : nids) {
    auto pids = mtx_helper.valve_loc_row2col[nid];
    for (Index pid : pids) {
      auto vid = std::make_pair<long, long>(nid, pid);
      auto valve = mtx_helper.id2vname[vid];
      valves.insert(valve);
      //        std::cout<<"n pick "<<valve<<std::endl;
    }
  }
  return valves;
}

void pipenetwork::ValveGraph::construct_iso_segs() {
  std::set<Index> unexplored_pids;
  for (const auto& index_link : mesh_links_->links()) {
    auto lid = index_link.first;
    unexplored_pids.insert(lid);
  }

  while (!unexplored_pids.empty()) {
    auto searching_pid = *unexplored_pids.begin();
    auto seg = get_single_seg(searching_pid);
    for (auto isolated_pid : seg.pids) {
      unexplored_pids.erase(isolated_pid);
      pid2segment_[isolated_pid] = seg;
    }
    iso_segs_.emplace_back(seg);
  }
}
