#ifndef PIPE_NETWORK_VALVE_GRAPH_COMPONENTS_H
#define PIPE_NETWORK_VALVE_GRAPH_COMPONENTS_H
#include "mesh_components.h"
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>

namespace pipenetwork {
namespace isolation {
//! Isolate Valve Property
//! name Valve name
//! on_node node name that the valve is closed to
//! on_pipe pipe name that the valve is on
struct ISOVProp {
  std::string name{"junction"};
  std::string on_node{"Null"};
  std::string on_pipe{"Null"};
};
struct IsoSeg {
  Index sid;
  std::set<Index> pids;
  std::set<Index> nids;
  std::set<Index> vids;
};
//! lookup tables for matrices used by the isolation finding algorithm
struct IsoMtxHelper {
  //! Index lookup table for the valve deficieny matrix (col2row)
  std::map<Index, std::vector<Index>> valve_def_col2row;
  //! Index lookup table for the valve deficieny matrix (row2col)
  std::map<Index, std::vector<Index>> valve_def_row2col;
  //! Index lookup table for the valve location matrix (col2row)
  std::map<Index, std::vector<Index>> valve_loc_col2row;
  //! Index lookup table for the valve location matrix (row2col)
  std::map<Index, std::vector<Index>> valve_loc_row2col;
};

struct IsoSegHelper {
  static const unsigned EIGEN_THRE{100};
  static constexpr double NON_ZERO_THRE{1E-8};
  static Eigen::SparseMatrix<double> shrink_mtx(
      Eigen::SparseMatrix<double>& matrix,
      const std::vector<pipenetwork::Index>& rowsToRemove,
      const std::vector<pipenetwork::Index>& colToRemove);

  static Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> small_matrix_eigen_info(
      Eigen::SparseMatrix<double>& matrix);

  static Spectra::SymEigsSolver<double, Spectra::SMALLEST_MAGN,
                                Spectra::SparseSymMatProd<double>>
      large_matrix_eigen_info(Eigen::SparseMatrix<double>& matrix);

  static Eigen::SparseMatrix<double> create_graph_laplacian(
      const Eigen::SparseMatrix<double>& adj_mtx);

  static std::vector<Index> find_zero_eval_loc(const Eigen::VectorXd& evals);

  static std::vector<Index> find_non_zero_evac_loc(const Eigen::VectorXd& evec);
};

class IsoValves {
 public:
  IsoValves() = default;

  IsoValves(const std::vector<ISOVProp>& iso_valve_props);

  //! Get the iso valves map
  const std::map<Index, ISOVProp>& iso_valves() { return iso_valves_; }

  //! Query the vid from  the position in the valve location matrix
  Index loc2vid(const std::pair<Index, Index> location) const {
    return loc2vid_.at(location);
  }

  //! Set location to vid map
  void loc2vid(std::pair<Index, Index> location, Index vid) {
    loc2vid_[location] = vid;
  }

  const Index nvalves() const { return iso_valves_.size(); }

 private:
  //! Internal valve id manager
  IndexManager vid_manager_;
  //! Valve properties map
  std::map<Index, ISOVProp> iso_valves_;
  //! valve name map with the position in the valve location matrix as the key
  std::map<std::pair<Index, Index>, Index> loc2vid_;
};

class IsoSegments {
 public:
  IsoSegments() = default;

  void construct_iso_segs(const IsoMtxHelper& mtx_helper,
                          pipenetwork::isolation::IsoValves& iso_valves,
                          std::set<Index>& pids);

  const isolation::IsoSeg get_iso_seg(Index pid);

  const tsl::ordered_map<Index, IsoSeg>& iso_segment() { return iso_segments_; }

  const Index nsegs() { return iso_segments_.size(); }

  Eigen::SparseMatrix<double>& seg_valve_mtx() { return seg_valve_mtx_; }

  Eigen::SparseMatrix<double>& seg_valve_adj_mtx() {
    return seg_valve_adj_mtx_;
  }

  void merge_segments(const std::vector<pipenetwork::Index>& broken_vids);

  std::vector<std::vector<Index>> get_segment_components(
      const std::vector<pipenetwork::Index>& segs2iso);

 private:
  //! Internal segment id manager
  IndexManager sid_manager_;
  //! Segment map
  tsl::ordered_map<Index, IsoSeg> iso_segments_;
  //! pid to sid
  tsl::ordered_map<Index, Index> pid2sid_;
  //! segment valve matrix
  Eigen::SparseMatrix<double> seg_valve_mtx_;
  //! Adjacency matrix for segments valve
  Eigen::SparseMatrix<double> seg_valve_adj_mtx_;
  //! merged segments. Key: from sid, value: to sid
  std::map<Index, Index> merged_sid_map_;

  void construct_seg_valve_mtx(pipenetwork::isolation::IsoValves& iso_valves);

  void construct_seg_valve_adj_mtx();

  //! get the corresponding isolation segment from a pipe id
  IsoSeg get_single_seg(const IsoMtxHelper& mtx_helper, Index pid);

  //! get the isolation valves that need to be closed from the corresponding
  //! segment pids and segment nids
  std::set<pipenetwork::Index> get_seg_valves(
      const pipenetwork::isolation::IsoMtxHelper& mtx_helper,
      const pipenetwork::isolation::IsoValves& iso_valves,
      const std::set<Index>& pids, const std::set<Index>& nids);

  void merge_two_segment(Index sid_from, Index sid_to);

  void update_seg_valve_mtx(const std::vector<pipenetwork::Index>& sids_remove,
                            const std::vector<pipenetwork::Index>& vids_remove);

  std::vector<Index> find_merging_sids(Index broken_vid);

  std::vector<Index> map2sid(std::vector<Index> eigen_ids);
};

}  // namespace isolation
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_VALVE_GRAPH_COMPONENTS_H
