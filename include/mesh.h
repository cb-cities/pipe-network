#ifndef PIPE_NETWORK_MESH_H
#define PIPE_NETWORK_MESH_H

#include "index_manager.h"
#include "junction.h"
#include "pipe.h"
#include "pump.h"
#include "reservoir.h"
#include "valve.h"

#include <array>
#include <exception>
//#include <cmath>
//#include <exception>
//#include <functional>
//#include <iostream>
//#include <memory>
#include <iostream>
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

  //! node name to id map
  tsl::ordered_map<std::string, Index> name2id;
  //! Junction map with respect to internal index
  tsl::ordered_map<Index, std::shared_ptr<Junction>> junctions;
  //! Reservoir map with respect to internal index
  tsl::ordered_map<Index, std::shared_ptr<Reservoir>> reservoirs;

 private:
  //! Internal id manager
  IndexManager nid_manager_;
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
  MeshLinks(std::vector<PipeProp>& pipe_props,
            std::vector<PumpProp>& pump_props,
            std::vector<ValveProp>& valve_props, const MeshNodes& mesh_nodes);

  //! Pipes map with respect to internal index
  tsl::ordered_map<Index, std::shared_ptr<Pipe>> pipes;
  //! Pumps map with respect to internal index
  tsl::ordered_map<Index, std::shared_ptr<Pump>> pumps;
  //! Valves map with respect to internal index
  tsl::ordered_map<Index, std::shared_ptr<Valve>> valves;

 private:
  //! Internal id manager
  IndexManager lid_manager_;

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

////! Mesh class
////! \brief Class for mesh that contains node and pipe pointers
class Mesh {

 public:
  //! Constructor with id
  //! \param[in] id mesh id
  explicit Mesh(std::string& name) : name_{name} {};

  //! Return name
  //! \retval name_ name of the mesh
  std::string name() const { return name_; }

  //  //! Create mesh from input object
  //  //! \param[in] IO pointer to the input object
  //  void create_mesh_from_inp(std::shared_ptr<Input>& IO);

  //! Create junction pointers
  //! \param[in] junc_props vector of junction properties
  void create_nodes(const std::vector<JunctionProp>& junc_props,
                    const std::vector<ReservoirProp>& res_props) {
    mesh_nodes_ = std::make_shared<MeshNodes>(junc_props, res_props);
  };

  //! Create Reservoir pointers
  //! \param[in] res_props vector of reservoir properties
  void create_links(std::vector<PipeProp>& pipe_props,
                    std::vector<PumpProp>& pump_props,
                    std::vector<ValveProp>& valve_props) {
    mesh_links_ =
        std::make_shared<MeshLinks>(pipe_props, pump_props, valve_props,*mesh_nodes_);
  };

  //! get nodes
  std::shared_ptr<MeshNodes> nodes() const { return mesh_nodes_; }
  //! get links
  std::shared_ptr<MeshLinks> links() const { return mesh_links_; }

  //! Print summary for the mesh
  //  void print_summary();

 private:
  //! mesh name
  std::string name_;
  //! mesh nodes
  std::shared_ptr<MeshNodes> mesh_nodes_;
  //! mesh links
  std::shared_ptr<MeshLinks> mesh_links_;
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_MESH_H
