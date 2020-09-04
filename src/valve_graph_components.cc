#include "valve_graph_components.h"

pipenetwork::isolation::IsoValves::IsoValves(
    const std::vector<ISOVProp>& iso_valve_props) {
  for (const auto& vprop : iso_valve_props) {
    Index vid = vid_manager_.create_index();
    iso_valves_[vid] = vprop;
    vname2vid_[vprop.name] = vid;
  }
}

const pipenetwork::isolation::IsoSeg
    pipenetwork::isolation::IsoSegments::pid2seg(pipenetwork::Index pid) {
  auto sid = pid2sid_[pid];
  return iso_segments_[sid];
}

void pipenetwork::isolation::IsoSegments::construct_iso_segs(
    const pipenetwork::isolation::IsoMtxHelper& mtx_helper,
    pipenetwork::isolation::IsoValves& iso_valves, std::set<Index>& pids) {
  while (!pids.empty()) {
    auto searching_pid = *pids.begin();
    auto sid = sid_manager_.create_index();
    auto seg = get_single_seg(mtx_helper, searching_pid);
    seg.vids = get_seg_valves(mtx_helper, iso_valves, seg.pids,
                              seg.nids);  // update vids
    for (auto isolated_pid : seg.pids) {
      pids.erase(isolated_pid);
      iso_segments_[sid] = seg;
    }
  }
  construct_seg_valve_mtx(iso_valves);
  construct_seg_valve_adj_mtx();
  update_pid2sid();
}

pipenetwork::isolation::IsoSeg
    pipenetwork::isolation::IsoSegments::get_single_seg(
        const pipenetwork::isolation::IsoMtxHelper& mtx_helper,
        pipenetwork::Index pid) {
  isolation::IsoSeg seg;
  seg.sid = sid_manager_.current_index();

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
    auto new_nids = mtx_helper.valve_def_col2row.at(pid);
    // row search (find new isolated pipes)
    while (!new_nids.empty()) {
      auto searching_nid = new_nids.back();
      new_nids.pop_back();
      seg.nids.insert(searching_nid);
      // adding new pids
      auto new_pids = mtx_helper.valve_def_row2col.at(searching_nid);
      for (auto new_pid : new_pids) {
        bool psearched = std::find(seg.pids.begin(), seg.pids.end(), new_pid) !=
                         seg.pids.end();
        if (!psearched) {
          pids.insert(new_pid);
        }
      }
    }
  }
  return seg;
}

std::set<pipenetwork::Index>
    pipenetwork::isolation::IsoSegments::get_seg_valves(
        const pipenetwork::isolation::IsoMtxHelper& mtx_helper,
        const pipenetwork::isolation::IsoValves& iso_valves,
        const std::set<Index>& pids, const std::set<Index>& nids) {
  std::set<Index> vids;
  // pick by pipes
  for (auto pid : pids) {
    if (mtx_helper.valve_loc_col2row.find(pid) !=
        mtx_helper.valve_loc_col2row.end()) {  // pipe may have no valve
      auto nids = mtx_helper.valve_loc_col2row.at(pid);
      for (auto nid : nids) {
        auto vloc = std::make_pair<long, long>(nid, pid);
        auto vid = iso_valves.loc2vid(vloc);
        vids.insert(vid);
      }
    }
  }
  // pick by nodes
  for (Index nid : nids) {
    if (mtx_helper.valve_loc_row2col.find(nid) !=
        mtx_helper.valve_loc_row2col.end()) {  // node may have no valve
      auto pids = mtx_helper.valve_loc_row2col.at(nid);
      for (Index pid : pids) {
        auto vloc = std::make_pair<long, long>(nid, pid);
        auto vid = iso_valves.loc2vid(vloc);
        vids.insert(vid);
      }
    }
  }
  return vids;
}

void pipenetwork::isolation::IsoSegments::construct_seg_valve_mtx(
    pipenetwork::isolation::IsoValves& iso_valves) {
  auto valves = iso_valves.iso_valves();
  std::vector<Eigen::Triplet<double>> graph_triplet;
  int on_val = 1;
  for (const auto& sid_seg : iso_segments_) {
    auto sid = sid_seg.first;
    auto seg = sid_seg.second;
    for (const auto& vid : seg.vids) {
      graph_triplet.emplace_back(sid, vid, on_val);
    }
  }
  seg_valve_mtx_.resize(iso_segments_.size(), iso_valves.nvalves());
  seg_valve_mtx_.setFromTriplets(graph_triplet.begin(), graph_triplet.end());
}

void pipenetwork::isolation::IsoSegments::construct_seg_valve_adj_mtx() {
  std::vector<Eigen::Triplet<double>> graph_triplet;
  auto ncols = seg_valve_mtx_.cols();
  for (int k = 0; k < ncols; ++k) {
    if (std::find(removed_vids_.begin(), removed_vids_.end(), k) ==
        removed_vids_.end()) {
      Eigen::VectorXd v_seg = seg_valve_mtx_.col(k);
      std::vector<Index> connected_segs;
      if (v_seg.sum() == 2) {
        for (Eigen::Index i = 0; i < v_seg.size(); ++i) {
          if (v_seg[i]) connected_segs.push_back(i);
        }
        graph_triplet.emplace_back(connected_segs[0], connected_segs[1], 1);
        graph_triplet.emplace_back(connected_segs[1], connected_segs[0], 1);
      }
    }
  }

  seg_valve_adj_mtx_.resize(iso_segments_.size(), iso_segments_.size());
  seg_valve_adj_mtx_.setFromTriplets(graph_triplet.begin(),
                                     graph_triplet.end());
}

void pipenetwork::isolation::IsoSegments::merge_segments(Index broken_vid) {
  std::vector<Index> sids_to_remove;
  auto sids = find_merging_sids(broken_vid);
  if (sids.size() == 2) {
    // one valve can only connect two segments
    auto sid_from = sids[1];
    auto sid_to = sids[0];
    merge_two_segment(sid_from, sid_to, broken_vid);
    sids_to_remove.emplace_back(sid_from);
  }
  update_seg_valve_mtx(sids_to_remove);
  reindex();  // reindex for consistency (one segment is removed)
  update_pid2sid();
}

std::vector<pipenetwork::Index>
    pipenetwork::isolation::IsoSegments::find_merging_sids(Index broken_vid) {
  Eigen::VectorXd v_seg = seg_valve_mtx_.col(broken_vid);
  std::vector<Index> segs_to_merge;
  for (Eigen::Index i = 0; i < v_seg.size(); ++i) {
    if (v_seg[i]) segs_to_merge.push_back(i);
  }
  return segs_to_merge;
}

void pipenetwork::isolation::IsoSegments::merge_two_segment(
    pipenetwork::Index sid_from, pipenetwork::Index sid_to,
    pipenetwork::Index vid) {
  auto& seg_from = iso_segments_[sid_from];
  auto& seg_to = iso_segments_[sid_to];

  seg_to.vids.insert(seg_from.vids.begin(), seg_from.vids.end());
  seg_to.pids.insert(seg_from.pids.begin(), seg_from.pids.end());
  seg_to.nids.insert(seg_from.nids.begin(), seg_from.nids.end());
  for (int j = 0; j < seg_valve_mtx_.cols(); j++) {
    seg_valve_mtx_.coeffRef(sid_to, j) += seg_valve_mtx_.coeffRef(sid_from, j);
  }

  // update states
  seg_to.vids.erase(vid);  // remove the broken valve
  removed_vids_.emplace_back(vid);
  iso_segments_.erase(sid_from);
}

void pipenetwork::isolation::IsoSegments::update_seg_valve_mtx(
    const std::vector<pipenetwork::Index>& sids_remove) {

  auto updated_seg_valve = pipenetwork::isolation::IsoSegHelper::shrink_mtx(
      seg_valve_mtx_, sids_remove, {});
  seg_valve_mtx_ = updated_seg_valve;
  //  std::cout << "Seg Valve mtx after!" << std::endl;
  //  std::cout << seg_valve_mtx_ << std::endl;
  construct_seg_valve_adj_mtx();
}

std::vector<std::vector<pipenetwork::isolation::IsoSeg>>
    pipenetwork::isolation::IsoSegments::get_segment_components(
        const std::set<pipenetwork::Index>& segs2iso) {

  // isolate segments by setting corresponding connections to 0
  auto n = seg_valve_adj_mtx_.cols();
  for (int i = 0; i < n; i++) {
    for (auto sid : segs2iso) {
      seg_valve_adj_mtx_.coeffRef(sid, i) = 0;
      seg_valve_adj_mtx_.coeffRef(i, sid) = 0;
    }
  }

  // construct graph laplacian
  Eigen::SparseMatrix<double> L =
      isolation::IsoSegHelper::create_graph_laplacian(seg_valve_adj_mtx_);
  // solve for eigen information
  Eigen::VectorXd eigen_vals;
  Eigen::MatrixXd eigen_vecs;
  auto eigensolver =
      pipenetwork::isolation::IsoSegHelper::small_matrix_eigen_info(L);
  eigen_vals = eigensolver.eigenvalues();
  eigen_vecs = eigensolver.eigenvectors();

  return einfo2components(eigen_vals, eigen_vecs);
}

std::vector<std::vector<pipenetwork::isolation::IsoSeg>>
    pipenetwork::isolation::IsoSegments::einfo2components(
        const Eigen::VectorXd& eigen_vals, const Eigen::MatrixXd& eigen_vecs) {
  auto valid_eids = isolation::IsoSegHelper::find_zero_eval_loc(eigen_vals);
  std::vector<std::vector<pipenetwork::isolation::IsoSeg>> components;
  for (auto valid_eid : valid_eids) {
    Eigen::VectorXd evec = eigen_vecs.col(valid_eid);
    std::vector<pipenetwork::isolation::IsoSeg> component;

    auto sids = isolation::IsoSegHelper::find_non_zero_evac_loc(evec);
    for (auto sid : sids) {
      component.emplace_back(iso_segments_[sid]);
    }
    components.emplace_back(component);
  }
  return components;
}

void pipenetwork::isolation::IsoSegments::reindex() {
  sid_manager_ = IndexManager();
  tsl::ordered_map<Index, IsoSeg> reindexed_segments;
  for (auto& sid_seg : iso_segments_) {
    auto seg = sid_seg.second;
    auto new_sid = sid_manager_.create_index();
    seg.sid = new_sid;
    reindexed_segments[new_sid] = seg;
  }
  iso_segments_ = reindexed_segments;
}

void pipenetwork::isolation::IsoSegments::update_pid2sid() {
  for (auto& sid_seg : iso_segments_) {
    auto seg = sid_seg.second;
    for (auto pid : seg.pids) {
      pid2sid_[pid] = seg.sid;
    }
  }
}

Eigen::SparseMatrix<double> pipenetwork::isolation::IsoSegHelper::shrink_mtx(
    Eigen::SparseMatrix<double>& matrix,
    const std::vector<pipenetwork::Index>& rowsToRemove,
    const std::vector<pipenetwork::Index>& colsToRemove) {

  Eigen::SparseMatrix<double> shrinked_mtx;
  shrinked_mtx.resize(matrix.rows() - rowsToRemove.size(),
                      matrix.cols() - colsToRemove.size());

  std::map<Index, Index> row_idx_map;
  for (int i = 0; i < matrix.innerSize(); ++i) {
    int decrement = 0;
    for (auto row2r : rowsToRemove) {
      if (i > row2r) {
        decrement -= 1;
      }
    }
    auto new_row_idx = i + decrement;
    row_idx_map[i] = new_row_idx;
  }

  Index col_new = 0;
  for (int k = 0; k < matrix.outerSize(); ++k) {
    if (std::find(colsToRemove.begin(), colsToRemove.end(), k) ==
        colsToRemove.end()) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, k); it; ++it) {
        auto row_old = it.row();  // row index
        auto val = it.value();
        if (std::find(rowsToRemove.begin(), rowsToRemove.end(), row_old) ==
                rowsToRemove.end() &&
            val != 0) {
          auto row_new = row_idx_map[row_old];
          shrinked_mtx.insert(row_new, col_new) = 1;
        }
      }
      col_new += 1;
    }
  }
  return shrinked_mtx;
}

Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>
    pipenetwork::isolation::IsoSegHelper::small_matrix_eigen_info(
        Eigen::SparseMatrix<double>& matrix) {

  auto dense_L = Eigen::MatrixXd(matrix);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(dense_L);
  if (eigensolver.info() != Eigen::Success) abort();
  return eigensolver;
}

Eigen::SparseMatrix<double>
    pipenetwork::isolation::IsoSegHelper::create_graph_laplacian(
        const Eigen::SparseMatrix<double>& adj_mtx) {
  std::vector<Eigen::Triplet<double>> graph_triplet;
  for (int i = 0; i < adj_mtx.innerSize(); i++) {
    auto ki = adj_mtx.row(i).sum();
    graph_triplet.emplace_back(i, i, ki);
  }
  Eigen::SparseMatrix<double> D;
  D.resize(adj_mtx.outerSize(), adj_mtx.outerSize());
  D.setFromTriplets(graph_triplet.begin(), graph_triplet.end());

  Eigen::SparseMatrix<double> L = D - adj_mtx;
  return L;
}

std::vector<pipenetwork::Index>
    pipenetwork::isolation::IsoSegHelper::find_zero_eval_loc(
        const Eigen::VectorXd& evals) {
  std::vector<pipenetwork::Index> valid_eids;
  for (int i = 0; i < evals.size(); i++) {
    auto val = evals[i];
    if (std::abs(val) < NON_ZERO_THRE) {
      valid_eids.emplace_back(i);
    }
  }
  return valid_eids;
}

std::vector<pipenetwork::Index>
    pipenetwork::isolation::IsoSegHelper::find_non_zero_evac_loc(
        const Eigen::VectorXd& evec) {

  std::vector<pipenetwork::Index> non_zero_ids;
  for (int i = 0; i < evec.size(); i++) {
    auto val = evec[i];
    if (std::abs(val) > NON_ZERO_THRE) {
      non_zero_ids.emplace_back(i);
    }
  }
  return non_zero_ids;
}
