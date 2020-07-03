
#include "valve_graph.h"

pipenetwork::ValveGraph::ValveGraph(
    const std::shared_ptr<pipenetwork::MeshNodes>& mesh_nodes,
    const std::shared_ptr<pipenetwork::MeshLinks>& mesh_links,
    const std::vector<pipenetwork::isolation::ISOVProp>& iso_valve_props)
    : MeshGraph(mesh_nodes, mesh_links), iso_valves_(iso_valve_props) {
  construct_node_pipe_mtx();
  construct_valve_loc_mtx();
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

void pipenetwork::ValveGraph::construct_iso_segs() {
  std::set<Index> unexplored_pids;
  for (const auto& index_link : mesh_links_->links()) {
    auto lid = index_link.first;
    unexplored_pids.insert(lid);
  }
  iso_segments_.construct_iso_segs(mtx_helper, iso_valves_, unexplored_pids);
}

// void pipenetwork::ValveGraph::find_segment_components(pipenetwork::Index sid)
// {
//
//  auto n = seg_valve_adj_mtx_.cols();
//  // isolate the segment to 0
//  for (int i = 0; i < n; i++) {
//    seg_valve_adj_mtx_.coeffRef(sid, i) = 0;
//    seg_valve_adj_mtx_.coeffRef(i, sid) = 0;
//  }
//
//  std::vector<Eigen::Triplet<double>> graph_triplet;
//  for (int i = 0; i < iso_segments_.nsegs(); i++) {
//    auto ki = seg_valve_adj_mtx_.row(i).sum();
//    graph_triplet.emplace_back(i, i, ki);
//  }
//  Eigen::SparseMatrix<double> D;
//  D.resize(iso_segments_.nsegs(), iso_segments_.nsegs());
//  D.setFromTriplets(graph_triplet.begin(), graph_triplet.end());
//
//  Eigen::SparseMatrix<double> L = D - seg_valve_adj_mtx_;
//  std::cout << L << std::endl;
//  if (n > pipenetwork::isolation::IsoSegHelper::EIGEN_THRE) {
//    auto eigensolver =
//        pipenetwork::isolation::IsoSegHelper::large_matrix_eigen_info(L);
//  } else {
//    auto eigensolver =
//        pipenetwork::isolation::IsoSegHelper::small_matrix_eigen_info(L);
//    std::cout << "The eigenvalues of A are:\n"
//              << eigensolver.eigenvalues() << std::endl;
//    std::cout << "Here's a matrix whose columns are eigenvectors of A \n"
//              << "corresponding to these eigenvalues:\n"
//              << eigensolver.eigenvectors() << std::endl;
//  }
//}
//
// Eigen::SparseMatrix<double> pipenetwork::ValveGraph::merge_segments(
//    Index vid, Eigen::SparseMatrix<double> seg_valve_mtx) {
//  Eigen::VectorXd v_seg = seg_valve_mtx.col(vid);
//  std::vector<Index> merge_pair;
//  for (Eigen::Index i = 0; i < v_seg.size(); ++i) {
//    if (v_seg[i]) merge_pair.push_back(i);
//  }
//  if (merge_pair.size() != 2) {
//    return seg_valve_mtx;
//  }
//
//  for (int j = 0; j < seg_valve_mtx.cols(); j++) {
//    seg_valve_mtx.coeffRef(merge_pair[0], j) +=
//        seg_valve_mtx.coeffRef(merge_pair[1], j);
//  }
//  auto merged_seg_valve_mtx =
//  pipenetwork::isolation::IsoSegHelper::shrink_mtx(
//      seg_valve_mtx, merge_pair[1], vid);
//  return merged_seg_valve_mtx;
//}
