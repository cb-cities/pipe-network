#ifndef PIPE_NETWORK_MESH_H_
#define PIPE_NETWORK_MESH_H_

#include <cmath>

#include <array>
#include <exception>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <vector>

#include "node.h"
#include "pipe.h"

//! Mesh class
//! \brief Class for mesh that contains node and pipe pointers
class Mesh {

 public:
  // Constructor with id
  //! \param[in] id mesh id
  Mesh(unsigned id) : id_{id} {};

  //! Destructor
  ~Mesh() = default;

  //! Return id
  //! \retval id_ id of the mesh
  unsigned id() const { return id_; }

  //! Create a nodal pointer based on input id and coordinates and check whether
  //! duplicate nodes exist \param[in] id vector of the node indices \param[in]
  //! coords vector of the coordinates of the nodes
  void create_node(std::vector<unsigned> id,
                   std::vector<Eigen::Vector3d> coords);

  //! Create a pipe pointer based on input idices of the pipe and the nodes at
  //! its ends and check whether duplicate pipe exist \param[in] pipeid index of
  //! the pipe \param[in] nodeid1 and nodeid2 indices of the nodes at pipe ends
  void create_pipe(std::vector<unsigned> pipeid,
                   std::vector<std::pair<unsigned, unsigned>> nodeids);

  //! Return the number of nodes in the mesh
  //! \retval nodes_.size() number of nodes
  unsigned nnodes() const { return nodes_.size(); }

  //! Return the number of pipes in the mesh
  //! \retval pipes_.size() number of pipes
  unsigned npipes() const { return pipes_.size(); }

  //! Check whether isolated node exists
  //! \retval the status to indicate whether isolated node exists
  bool isolated_node();

  //! Return coordinates of all the nodes in the mesh
  //! \retval nodal_coordinates coordinates of all the nodes
  std::vector<Eigen::Vector3d> nodal_coordinates();

 private:
  //! mesh id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! nodal id and corresponding nodal pointer
  std::map<unsigned, std::shared_ptr<pipenetwork::Node>> nodes_;
  //! pipe id and corresponding pipe pointer
  std::map<unsigned, std::unique_ptr<pipenetwork::Pipe>> pipes_;
};

#endif  // PIPE_NETWORK_MESH_H_
