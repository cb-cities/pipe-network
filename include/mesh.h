#ifndef PIPE_NETWORK_MESH_H
#define PIPE_NETWORK_MESH_H

#include "pipe.h"
#include "pump.h"
#include "junction.h"
#include "reservoir.h"

#include <array>
//#include <cmath>
//#include <exception>
//#include <functional>
//#include <iostream>
//#include <memory>
#include <tsl/ordered_map.h>
#include <tuple>
#include <vector>


namespace pipenetwork {

//! Mesh class
//! \brief Class for mesh that contains node and pipe pointers
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
  void create_junctions(const std::vector<JunctionProp>& junc_props);

  //! Create Reservoir pointers
  //! \param[in] res_props vector of reservoir properties
  void create_reservoirs(const std::vector<ReservoirProp>& res_props);

  //! Create Pipe pointers
  //! \param[in]  pipe_props vector of pipe properties
  void create_pipes(std::vector<PipeProp>& pipe_props);

  //! Create Pump pointers
  //! \param[in]  pump_props vector of pump properties
  void create_pumps(std::vector<PumpProp>& pump_props);

  //! Create Valve pointers
  //! \param[in]  valve_props vector of valve properties
  void create_valve(std::vector<Valve_prop>& valve_props);

  //! get all nodes map
  tsl::ordered_map<std::string, std::shared_ptr<pipenetwork::Node>> nodes()
      const {
    return nodes_;
  }
  //! get links map
  std::vector<std::shared_ptr<pipenetwork::Link>> links() const {
    return links_;
  }

  //! Iterate over nodes
  //! \tparam Toper Callable object typically a baseclass functor
  template <typename Toper>
  void iterate_over_nodes(Toper oper) {
    std::for_each(
        nodes_.cbegin(), nodes_.cend(),
        [=](std::pair<std::string, std::shared_ptr<pipenetwork::Node>> node) {
          oper(node.second);
        });
  };

  //! Iterate over links
  //! \tparam Toper Callable object typically a baseclass functor
  template <typename Toper>
  void iterate_over_links(Toper oper) {
    std::for_each(links_.cbegin(), links_.cend(),
                  [=](std::shared_ptr<pipenetwork::Link> link) { oper(link); });
  };

  //! Print summary for the mesh
  void print_summary();

  //! Return number of nodes in the network
  //! \retval nnode_ number of nodes in the network
  unsigned nnodes() const { return nodes_.size(); }

  //! Return number of links in the network
  //! \retval nnode_ number of pipes in the network
  unsigned nlinks() const { return links_.size(); }

  //! Return number of pipes in the network
  //! \retval nnode_ number of pipes in the network
  unsigned npipes() const { return npipes_; }

  //! Return number of pumps in the network
  //! \retval nnode_ number of pumps in the network
  unsigned npumps() const { return npumps_; }

  //! Return number of valves in the network
  //! \retval nnode_ number of valves in the network
  unsigned nvalves() const { return nvalves_; }

  //! Return number of junctions in the network
  //! \retval nnode_ number of valves in the network
  unsigned njunctions() const { return njunctions_; }

  //! Return number of sources in the network
  //! \retval nnode_ number of valves in the network
  unsigned nsources() const { return nsrcs_; }

 private:
  //! mesh name
  std::string name_;

  //! number of pumps
  unsigned npumps_{0};
  //! number of pipes
  unsigned npipes_{0};
  //! number of valves
  unsigned nvalves_{0};
  //! number of junctions
  unsigned njunctions_{0};
  //! number of sources
  unsigned nsrcs_{0};
  //! nodal id and corresponding nodal pointer
  tsl::ordered_map<std::string, std::shared_ptr<pipenetwork::Node>> nodes_;
  //! vector of links
  std::vector<std::shared_ptr<pipenetwork::Link>> links_;
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_MESH_H
