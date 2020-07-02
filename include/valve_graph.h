#ifndef PIPE_NETWORK_VALVE_GRAPH_H
#define PIPE_NETWORK_VALVE_GRAPH_H

#include "mesh_components.h"
#include "valve_graph_components.h"

#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>

namespace pipenetwork {

class ValveGraph : MeshGraph {
 public:
  ValveGraph(const std::shared_ptr<MeshNodes>& mesh_nodes,
             const std::shared_ptr<MeshLinks>& mesh_links,
             const std::vector<isolation::ISOVProp>& iso_valve_props);

  Eigen::SparseMatrix<int>& node_pipe_mtx() { return node_pipe_mtx_; }

  Eigen::SparseMatrix<int>& valve_loc_mtx() { return valve_loc_mtx_; }

  Eigen::SparseMatrix<double>& seg_valve_mtx() { return seg_valve_mtx_; }

  Eigen::SparseMatrix<double>& seg_valve_adj_mtx() {
    return seg_valve_adj_mtx_;
  }

  void find_segment_components(Index sid);

  //! get the corresponding isolation segment from a pipe
  const isolation::IsoSeg get_iso_seg(Index pid) {
    return iso_segments_.get_iso_seg(pid);
  }

  //! get number of segments
  const Index nsegs() { return iso_segments_.nsegs(); }

  Eigen::SparseMatrix<double> merge_segments(Index vid) {
    return merge_segments(vid, seg_valve_mtx_);
  }

 private:
  //! Isolation valves in the valve graph
  isolation::IsoValves iso_valves_;
  //! Isolation segments in the valve graph
  isolation::IsoSegments iso_segments_;

  //! node pipe matrix
  Eigen::SparseMatrix<int> node_pipe_mtx_;
  //! valve location matrix
  Eigen::SparseMatrix<int> valve_loc_mtx_;
  //! valve deficiency matrix
  Eigen::SparseMatrix<int> valve_def_mtx_;
  //! segment valve matrix
  Eigen::SparseMatrix<double> seg_valve_mtx_;
  //! Adjacency matrix for segments valve
  Eigen::SparseMatrix<double> seg_valve_adj_mtx_;
  //! lookup tables for the stored matrices
  isolation::IsoMtxHelper mtx_helper;

  //! construct the node_pipe_mtx_
  void construct_node_pipe_mtx();

  //! construct the valve location matrix
  void construct_valve_loc_mtx();

  //! construct the valve defeciency matrix
  void construct_valve_def_mtx();

  //! construct the col2row, row2col look up tables
  void construct_vdef_idx_table();

  //! construct isolation segments for the network
  void construct_iso_segs();

  void construct_seg_valve_mtx();

  void construct_seg_valve_adj_mtx();

  Eigen::SparseMatrix<double> merge_segments(
      Index vid, Eigen::SparseMatrix<double> seg_valve_mtx);
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_VALVE_GRAPH_H
