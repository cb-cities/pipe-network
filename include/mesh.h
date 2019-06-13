#ifndef PIPE_NETWORK_MESH_H
#define PIPE_NETWORK_MESH_H

#include <array>
#include <cmath>
#include <exception>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <tuple>
#include <vector>

// TBB
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>

#include "junction.h"
#include "pipe.h"
#include "reservoir.h"

namespace pipenetwork {

//! Mesh class
//! \brief Class for mesh that contains node and pipe pointers
class Mesh {

 public:
  //! Constructor with id
  //! \param[in] id mesh id
  explicit Mesh(unsigned id) : id_{id} {};

  //! Destructor
  ~Mesh() = default;

  //! Return id
  //! \retval id_ id of the mesh
  unsigned id() const { return id_; }

  //! Create junction pointers
  //! \param[in] ids vector of junction id
  //! \param[in] elevations vector of elevation for the junction
  //! \param[in] demands base vector of demand for the junction
  //! \param[in] leak_diameters vector of diameter of the leak hole for the
  //! junction
  void create_junctions(const std::vector<Index>& ids,
                        const std::vector<double>& elevations,
                        const std::vector<double>& demands,
                        const std::vector<double>& leak_diameters);

  //! Create Reservoir pointers
  //! \param[in] ids reservoir id
  //! \param[in] heads base head for reservoirs
  void create_reservoirs(const std::vector<Index>& ids,
                         const std::vector<double>& heads);

  //! Create Pipe pointers
  //! \param[in] nodeids pair of end node ids for the pipe
  //! \param[in] length, length of the pipe
  //! \param[in] diameter, diameter of the pipe
  //! \param[in] roughness, roughness of the pipe
  //! \param[in] status, status of the pipe (open or close)
  bool create_pipes(const std::vector<Index>& ids,
                    const std::vector<std::pair<Index, Index>>& nodeids,
                    const std::vector<double>& length,
                    const std::vector<double>& diameter,
                    const std::vector<double>& roughness,
                    const std::vector<Pipe_status>& status);

  //! get all nodes map
  std::map<Index, std::shared_ptr<pipenetwork::Node>> nodes() const {
    return nodes_;
  }
  //! get connected nodes map
  std::map<Index, std::shared_ptr<pipenetwork::Node>> connect_nodes() const {
    return connected_nodes_;
  }
  //! get links map
  std::map<Index, std::shared_ptr<pipenetwork::Link>> links() const {
    return links_;
  }

  //! Iterate over nodes
  //! \tparam Toper Callable object typically a baseclass functor
  template <typename Toper>
  void iterate_over_nodes(Toper oper) {
    std::for_each(
        connected_nodes_.cbegin(), connected_nodes_.cend(),
        [=](std::pair<Index, std::shared_ptr<pipenetwork::Node>> node) {
          oper(node.second);
        });
  };

  //! Iterate over links
  //! \tparam Toper Callable object typically a baseclass functor
  template <typename Toper>
  void iterate_over_links(Toper oper) {
    std::for_each(
        links_.cbegin(), links_.cend(),
        [=](std::pair<Index, std::shared_ptr<pipenetwork::Link>> link) {
          oper(link.second);
        });
  };

  //! Print summary for the mesh
  void print_summary();

  //! Return number of nodes in the network
  //! \retval nnode_ number of nodes in the network
  unsigned nnodes() const { return connected_nodes_.size(); }

  //! Return number of pipes in the network
  //! \retval nnode_ number of pipes in the network
  unsigned nlinks() const { return links_.size(); }

 private:
  //! mesh id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! nodal id and corresponding nodal pointer
  std::map<Index, std::shared_ptr<pipenetwork::Node>> nodes_;
  //! pipe id and corresponding pipe pointer
  std::map<Index, std::shared_ptr<pipenetwork::Link>> links_;
  //! nodal id and corresponding nodal pointer
  std::map<Index, std::shared_ptr<pipenetwork::Node>> connected_nodes_;
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_MESH_H
