#ifndef PIPE_NETWORK_MESH_H_
#define PIPE_NETWORK_MESH_H_

#include <cmath>

#include <array>
#include <exception>
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
  ~Mesh() {
    nodes_.clear();
    pipes_.clear();
  };

  //! Return id
  //! \retval id_ id of the mesh
  unsigned id() const { return id_; }

  //! Add a node pointer
  //! \param[in] nodeptr pointer to a node
  void add_node(std::shared_ptr<Node>& node) { nodes_.emplace_back(node); }

  //! Return a node pointer for a given index
  //! \param[in] id index of the desired node
  //! \retval nodes_.at(id) node pointer at the given index
  // std::shared_ptr<Node> nodes(unsigned id) const { return nodes_.at(id); }

  //! Add a pipe pointer
  //! \param[in] pipeptr pointer to a pipe
  void add_pipe(std::shared_ptr<Pipe>& pipe) {
    unsigned count = 0;
    for (const auto& node : nodes_)
      if (pipe->isnode(node)) count++;
    if (count == 2)
      pipes_.emplace_back(pipe);
    else if (count == 1)
      throw std::runtime_error("one node doesn't exist, cannot add the pipe");
    else
      throw std::runtime_error("nodes don't exist, cannot add the pipe");
  }

  //! Return a pipe pointer for a given index
  //! \param[in] id index of the desired pipe
  //! \retval pipes_.at(id) pipe pointer at the given index
  // std::shared_ptr<Pipe> pipes(unsigned id) const { return pipes_.at(id); }

  //! Return the number of nodes in the mesh
  //! \retval nodes_.size() number of nodes
  unsigned nnodes() const { return nodes_.size(); }

  //! Return the number of pipes in the mesh
  //! \retval pipes_.size() number of pipes
  unsigned npipes() const { return pipes_.size(); }

  //! Check whether redundant node exists
  void redundant_node() {
    bool connect = false;
    for (const auto& node : nodes_) {
      for (const auto& pipe : pipes_)
        if (pipe->isnode(node)) connect = true;
      if (connect == false) {
        throw std::runtime_error("redundant node exist, check input");
        break;
      }
    }
  }

  //! Return coordinates of all the nodes in the mesh
  //! \retval nodal_coordinates coordinates of all the nodes
  std::vector<Eigen::Vector3d> nodal_coordinates() {
    std::vector<Eigen::Vector3d> nodal_coordinates;
    for (const auto& node : nodes_)
      nodal_coordinates.emplace_back(node->coordinates());
    return nodal_coordinates;
  }

  //! Return degree centrality (number of pipe connected to the node) of a given
  //! node param[in] id index of the interested node \retval degree_centrality
  //! degree centrality of the node
  unsigned degree_centrality(unsigned id) {
    unsigned degree_centrality = 0;
    for (const auto& pipe : pipes_)
      if (pipe->isnode(nodes_.at(id))) degree_centrality++;
    return degree_centrality;
  }

 private:
  //! mesh id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! node pointers
  std::vector<std::shared_ptr<Node>> nodes_;
  //! pipe pointers
  std::vector<std::shared_ptr<Pipe>> pipes_;
};

#endif  // PIPE_NETWORK_MESH_H_
