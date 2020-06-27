#ifndef PIPE_NETWORK_VALVE_GRAPH_H
#define PIPE_NETWORK_VALVE_GRAPH_H

#include "mesh_components.h"
namespace pipenetwork {

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
  std::set<Index> pids;
  std::set<Index> nids;
  std::set<std::string> vnames;
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

  //! valve name map with the position in the valve location matrix as the key
  std::map<std::pair<Index, Index>, std::string> id2vname;
};

class ValveGraph : MeshGraph {
 public:
  ValveGraph(const std::shared_ptr<MeshNodes>& mesh_nodes,
             const std::shared_ptr<MeshLinks>& mesh_links,
             const std::vector<ISOVProp>& iso_valve_props);

  Eigen::SparseMatrix<int>& node_pipe_mtx() { return node_pipe_mtx_; }

  Eigen::SparseMatrix<int>& valve_loc_mtx() { return valve_loc_mtx_; }

  //! get the corresponding isolation segment from a pipe
  const IsoSeg get_iso_seg(Index pid) { return pid2segment_.at(pid); }

  //! get all isolation segments of the network
  const std::vector<IsoSeg> get_iso_segs() { return iso_segs_; }

 private:
  //! node pipe matrix
  Eigen::SparseMatrix<int> node_pipe_mtx_;
  //! valve location matrix
  Eigen::SparseMatrix<int> valve_loc_mtx_;
  //! valve deficiency matrix
  Eigen::SparseMatrix<int> valve_def_mtx_;
  //! lookup tables for the stored matrices
  IsoMtxHelper mtx_helper;

  std::vector<IsoSeg> iso_segs_;
  std::map<Index, IsoSeg> pid2segment_;

  //! construct the node_pipe_mtx_
  void construct_node_pipe_mtx();

  //! construct the valve location matrix
  void construct_valve_loc_mtx(const std::vector<ISOVProp>& iso_valve_props);

  //! construct the valve defeciency matrix
  void construct_valve_def_mtx();

  //! construct the col2row, row2col look up tables
  void construct_vdef_idx_table();

  //! get the isolation valves that need to be closed from corresponding pids
  //! and nids
  std::set<std::string> get_seg_valves(std::set<Index>& pids,
                                       std::set<Index>& nids);

  //! get the corresponding isolation segment from a pipe id
  IsoSeg get_single_seg(Index pid);

  //! construct isolation segments for the network
  void construct_iso_segs();
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_VALVE_GRAPH_H
