#include "valve_graph_components.h"

pipenetwork::isolation::IsoValves::IsoValves(
    const std::vector<ISOVProp>& iso_valve_props) {
  for (const auto& vprop : iso_valve_props) {
    Index vid = vid_manager_.create_index();
    iso_valves_[vid] = vprop;
  }
}

const pipenetwork::isolation::IsoSeg
    pipenetwork::isolation::IsoSegments::get_iso_seg(pipenetwork::Index pid) {
  auto sid = pid2sid_[pid];

  if (iso_segments_.find(sid) == iso_segments_.end()) {
    sid = merged_sid_map_[sid];  // segement has been merged, redirect to the
                                 // merged one
  }

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
      pid2sid_[isolated_pid] = sid;
    }
  }
  construct_seg_valve_mtx(iso_valves);
  construct_seg_valve_adj_mtx();
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
    auto nids = mtx_helper.valve_loc_col2row.at(pid);
    for (auto nid : nids) {
      auto vloc = std::make_pair<long, long>(nid, pid);
      auto vid = iso_valves.loc2vid(vloc);
      vids.insert(vid);
    }
  }
  // pick by nodes
  for (Index nid : nids) {
    auto pids = mtx_helper.valve_loc_row2col.at(nid);
    for (Index pid : pids) {
      auto vloc = std::make_pair<long, long>(nid, pid);
      auto vid = iso_valves.loc2vid(vloc);
      vids.insert(vid);
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
  seg_valve_adj_mtx_.resize(iso_segments_.size(), iso_segments_.size());
  seg_valve_adj_mtx_.setFromTriplets(graph_triplet.begin(),
                                     graph_triplet.end());
}

void pipenetwork::isolation::IsoSegments::merge_segments(
    const std::vector<pipenetwork::Index>& broken_vids) {
  std::vector<Index> sids_to_remove;
  for (const auto broken_vid : broken_vids) {
    auto sids = find_merging_sids(broken_vid);
    if (sids.size() == 2) {
      // one valve can only connect two segments
      auto sid_from = sids[1];
      auto sid_to = sids[0];
      merge_two_segment(sid_from, sid_to);
      sids_to_remove.emplace_back(sid_from);
    }
  }
  update_seg_valve_mtx(sids_to_remove, broken_vids);
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
    pipenetwork::Index sid_from, pipenetwork::Index sid_to) {
  auto& seg_from = iso_segments_[sid_from];
  auto& seg_to = iso_segments_[sid_to];

  seg_to.vids.insert(seg_from.vids.begin(), seg_from.vids.end());
  seg_to.pids.insert(seg_from.pids.begin(), seg_from.pids.end());
  seg_to.nids.insert(seg_from.nids.begin(), seg_from.nids.end());
  for (int j = 0; j < seg_valve_mtx_.cols(); j++) {
    seg_valve_mtx_.coeffRef(sid_to, j) += seg_valve_mtx_.coeffRef(sid_from, j);
  }

  // update states
  iso_segments_.erase(sid_from);
  merged_sid_map_[sid_from] = sid_to;
}

void pipenetwork::isolation::IsoSegments::update_seg_valve_mtx(
    const std::vector<pipenetwork::Index>& sids_remove,
    const std::vector<pipenetwork::Index>& vids_remove) {

  auto updated_seg_valve = pipenetwork::isolation::IsoSegHelper::shrink_mtx(
      seg_valve_mtx_, sids_remove, vids_remove);
  seg_valve_mtx_ = updated_seg_valve;
}

std::vector<std::vector<pipenetwork::Index>>
    pipenetwork::isolation::IsoSegments::get_segment_components(
        const std::vector<pipenetwork::Index>& segs2iso) {

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

  Eigen::VectorXd eigen_vals;
  Eigen::MatrixXd eigen_vecs;
  if (n > pipenetwork::isolation::IsoSegHelper::EIGEN_THRE) {
    auto eigensolver =
        pipenetwork::isolation::IsoSegHelper::large_matrix_eigen_info(L);
    eigen_vals = eigensolver.eigenvalues();
    eigen_vecs = eigensolver.eigenvectors();
  } else {
    auto eigensolver =
        pipenetwork::isolation::IsoSegHelper::small_matrix_eigen_info(L);
    eigen_vals = eigensolver.eigenvalues();
    eigen_vecs = eigensolver.eigenvectors();
  }

  auto valid_eids = isolation::IsoSegHelper::find_zero_eval_loc(eigen_vals);

  std::vector<std::vector<pipenetwork::Index>> components;
  for (auto valid_eid : valid_eids) {
    Eigen::VectorXd evec = eigen_vecs.col(valid_eid);
    auto eigen_ids = isolation::IsoSegHelper::find_non_zero_evac_loc(evec);
    components.emplace_back(map2sid(eigen_ids));
  }
  return components;
}

std::vector<pipenetwork::Index> pipenetwork::isolation::IsoSegments::map2sid(
    std::vector<pipenetwork::Index> eigen_ids) {
  long i = 0;
  std::vector<pipenetwork::Index> sids;
  for (const auto& sid_seg : iso_segments_) {
    if (std::find(eigen_ids.begin(), eigen_ids.end(), i) != eigen_ids.end()) {
      sids.emplace_back(sid_seg.first);
    }
    i++;
  }
  return sids;
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

Spectra::SymEigsSolver<double, Spectra::SMALLEST_MAGN,
                       Spectra::SparseSymMatProd<double>>
    pipenetwork::isolation::IsoSegHelper::large_matrix_eigen_info(
        Eigen::SparseMatrix<double>& matrix) {

  Spectra::SparseSymMatProd<double> op(matrix);
  Spectra::SymEigsSolver<double, Spectra::SMALLEST_MAGN,
                         Spectra::SparseSymMatProd<double>>
      eigs(&op, EIGEN_THRE, 6);
  eigs.init();
  int nconv = eigs.compute();
  if (eigs.info() != Spectra::SUCCESSFUL) abort();
  return eigs;
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
    if (val < NON_ZERO_THRE) {
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
    if (val > NON_ZERO_THRE) {
      non_zero_ids.emplace_back(i);
    }
  }
  return non_zero_ids;
}
