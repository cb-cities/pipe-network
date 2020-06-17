#ifndef PIPE_NETWORK_MESH_COMPONENTS_H
#define PIPE_NETWORK_MESH_COMPONENTS_H
#include "index_manager.h"
#include "junction.h"
#include "pipe.h"
#include "pump.h"
#include "reservoir.h"
#include "valve.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>
#include <exception>
#include <iostream>
#include <set>
#include <tsl/ordered_map.h>
#include <tuple>
#include <vector>

namespace pipenetwork {
//! MeshNodes class
//! \brief Nodes information for a mesh
class MeshNodes {
 public:
  MeshNodes() = default;
  //! Constructor for different mesh nodes
  //! \param[in] junc_props junction properties
  //! \param[in] res_props reservoir properties
  MeshNodes(const std::vector<JunctionProp>& junc_props,
            const std::vector<ReservoirProp>& res_props);

  //! Get node ptr with node name
  //! \param[in] node_name node name
  //! \return node node ptr
  std::shared_ptr<Node> get_node(const std::string& node_name) const;

  //! Get junctions map
  tsl::ordered_map<Index, std::shared_ptr<Junction>>& junctions() {
    return junctions_;
  }

  //! Get reservoirs map
  tsl::ordered_map<Index, std::shared_ptr<Reservoir>>& reservoirs() {
    return reservoirs_;
  }
  //! Get nodes map
  const tsl::ordered_map<Index, std::shared_ptr<Node>>& nodes() {
    return nodes_;
  }

  //! number of junctions
  Index njunctions() const { return junctions_.size(); }

  //! number of reservoirs
  Index nreservoirs() const { return reservoirs_.size(); }

  //! number of nodes
  Index nnodes() const { return nodes_.size(); }

 private:
  //! Internal id manager
  IndexManager nid_manager_;
  //! node name to id map
  tsl::ordered_map<std::string, Index> name2id_;
  //! Junction map with respect to internal index
  tsl::ordered_map<Index, std::shared_ptr<Junction>> junctions_;
  //! Reservoir map with respect to internal index
  tsl::ordered_map<Index, std::shared_ptr<Reservoir>> reservoirs_;
  //! Nodes map with respect to internal index
  tsl::ordered_map<Index, std::shared_ptr<Node>> nodes_;

  //! Add list of nodes
  //! \param[in] props vector of node properties
  template <typename Prop>
  void add_nodes(const std::vector<Prop>& props);

  //! Add a single node (junction)
  //! \param[in] junc_prop junction property
  void add_node(const JunctionProp& junc_prop);

  //! Add a single node (reservoir)
  //! \param[in] res_prop reservoir property
  void add_node(const ReservoirProp& res_prop);
};

//! Meshlink class
//! \brief links information for a mesh
class MeshLinks {
 public:
  MeshLinks() = default;
  //! Constructor for different mesh links
  //! \param[in] pipe_props pipe properties
  //! \param[in] pump_props pump properties
  //! \param[in] valve_props valve properties
  //! \param[in] mesh_nodes information of nodes inside the mesh
  MeshLinks(const std::vector<PipeProp>& pipe_props,
            const std::vector<PumpProp>& pump_props,
            const std::vector<ValveProp>& valve_props,
            const MeshNodes& mesh_nodes);

  //! Get pipes map
  tsl::ordered_map<Index, std::shared_ptr<Pipe>>& pipes() { return pipes_; }

  //! Get pumps map
  tsl::ordered_map<Index, std::shared_ptr<Pump>>& pumps() { return pumps_; }

  //! Get valves map
  tsl::ordered_map<Index, std::shared_ptr<Valve>>& valves() { return valves_; }

  //! Get links map
  const tsl::ordered_map<Index, std::shared_ptr<Link>>& links() {
    return links_;
  }

  //! number of pipes
  Index npipes() const { return pipes_.size(); }

  //! number of pumps
  Index npumps() const { return pumps_.size(); }

  //! number of valves
  Index nvalves() const { return valves_.size(); }

  //! number of nodes
  Index nlinks() const { return links_.size(); }

 private:
  //! Internal id manager
  IndexManager lid_manager_;
  //! Pipes map with respect to internal index
  tsl::ordered_map<Index, std::shared_ptr<Pipe>> pipes_;
  //! Pumps map with respect to internal index
  tsl::ordered_map<Index, std::shared_ptr<Pump>> pumps_;
  //! Valves map with respect to internal index
  tsl::ordered_map<Index, std::shared_ptr<Valve>> valves_;
  //! Link map with respect to internal index
  tsl::ordered_map<Index, std::shared_ptr<Link>> links_;
  //! Add list of links
  //! \param[in] props vector of link properties
  //! \param[in] mesh_nodes information of nodes inside the mesh
  template <typename Prop>
  void add_links(const std::vector<Prop>& props, const MeshNodes& mesh_nodes);

  //! Add a single link (pipe)
  //! \param[in] pipe_prop pipe property
  void add_link(const std::shared_ptr<Node>& node1,
                const std::shared_ptr<Node>& node2, const PipeProp& pipe_prop);

  //! Add a single link (pump)
  //! \param[in] pipe_prop pump property
  void add_link(const std::shared_ptr<Node>& node1,
                const std::shared_ptr<Node>& node2, const PumpProp& pump_prop);

  //! Add a single link (valve)
  //! \param[in] valve_prop valve property
  void add_link(const std::shared_ptr<Node>& node1,
                const std::shared_ptr<Node>& node2,
                const ValveProp& valve_prop);
};

//! Mesh graph class
//! \brief graph information for a mesh
class MeshGraph {
 public:
  MeshGraph(const std::shared_ptr<MeshNodes>& mesh_nodes,
            const std::shared_ptr<MeshLinks>& mesh_links)
      : mesh_nodes_{mesh_nodes}, mesh_links_{mesh_links} {
    compute_graph_info_();
  };

  //! Get adjacency matrix of the mesh
  const Eigen::SparseMatrix<int>& adjacency_matrix() const { return A_; }
  //! Get node to link map
  const std::map<Index, std::vector<Index>>& node2link_map() const {
    return node2link_;
  }
  //! Get node degrees
  std::vector<unsigned> ndegree() const { return ndegree_; }
  //! BFS
  Eigen::VectorXd bfs(Index nid);

 private:
  //! mesh nodes
  std::shared_ptr<MeshNodes> mesh_nodes_;
  //! mesh links
  std::shared_ptr<MeshLinks> mesh_links_;
  //! adjacency matrix
  Eigen::SparseMatrix<int> A_;
  //! node id to link id
  std::map<Index, std::vector<Index>> node2link_;
  //! node degrees
  std::vector<unsigned> ndegree_;

  //! compute adjacency matrix
  void compute_graph_info_();

  //! compute node connectivity degrees
  void compute_node_degrees_();
};
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_MESH_COMPONENTS_H
