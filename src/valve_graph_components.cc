#include <valve_graph_components.h>

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
  return iso_segments_[sid];
}

void pipenetwork::isolation::IsoSegments::construct_iso_segs(
    const pipenetwork::isolation::IsoMtxHelper& mtx_helper,
    const pipenetwork::isolation::IsoValves& iso_valves,
    std::set<Index>& pids) {
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
