#ifndef PIPE_NETWORK_VALVE_GRAPH_COMPONENTS_H
#define PIPE_NETWORK_VALVE_GRAPH_COMPONENTS_H
#include "mesh_components.h"

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

  const Index nvalves() { return iso_valves_.size(); }

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
                          const pipenetwork::isolation::IsoValves& iso_valves,
                          std::set<Index>& pids);

  const isolation::IsoSeg get_iso_seg(Index pid);

  const std::map<Index, IsoSeg>& iso_segment() { return iso_segments_; }

  const Index nsegs() { return iso_segments_.size(); }

 private:
  //! Internal segment id manager
  IndexManager sid_manager_;
  //! Segment map
  std::map<Index, IsoSeg> iso_segments_;
  //! pid to sid
  std::map<Index, Index> pid2sid_;

  //! get the corresponding isolation segment from a pipe id
  IsoSeg get_single_seg(const IsoMtxHelper& mtx_helper, Index pid);

  //! get the isolation valves that need to be closed from corresponding pids
  //! and nids
  std::set<pipenetwork::Index> get_seg_valves(
      const pipenetwork::isolation::IsoMtxHelper& mtx_helper,
      const pipenetwork::isolation::IsoValves& iso_valves,
      const std::set<Index>& pids, const std::set<Index>& nids);
};

}  // namespace isolation
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_VALVE_GRAPH_COMPONENTS_H
