#ifndef PIPE_NETWORK_MESH_H_
#define PIPE_NETWORK_MESH_H_

#include <cmath>

#include <array>
#include <exception>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include "node.h"
#include "pipe.h"
#include "settings.h"

namespace pipenetwork {

//! Mesh class
//! \brief Class for mesh that contains node and pipe pointers
class Mesh {

 public:
  // Constructor with id
  //! \param[in] id mesh id
  explicit Mesh(unsigned id) : id_{id} {}  ;

  //! Destructor
  ~Mesh() = default;

  //! Return id
  //! \retval id_ id of the mesh
  unsigned id() const { return id_; }

  //! Create nodal pointers and assign indices based on coordinates
  //! \param[in] coords vector of the coordinates of the nodes
  void create_nodes(const std::vector<Eigen::Vector3d>& coords);

  //! Create a pipe pointers and assign indices based on the nodes at its ends
  //! \param[in] nodeid1 and nodeid2 indices of the nodes at pipe ends
  //! \retval status to check whether all input pipe created successfully
  bool create_pipes(const std::vector<std::pair<Index, Index>>& nodeids,
                    const std::vector<double>& diameter,
                    const std::vector<double>& roughness,
                    const std::vector<bool>& pipe_status);

  //! Return the number of nodes in the mesh
  //! \retval nodes_.size() number of nodes
  unsigned nnodes() const { return nodes_.size(); }

  //! Return the number of pipes in the mesh
  //! \retval pipes_.size() number of pipes
  unsigned npipes() const { return pipes_.size(); }

  //! Remove unconnected nodes from the mesh
  void remove_unconnected_nodes();

  //! Initialize discharges in pipes
  void initialize_pipe_discharge(
      const std::vector<std::pair<Index, double>>& init_discharge =
          std::vector<std::pair<Index, double>>());
  void initialize_pipe_discharge(double init_discharge);

  //! Assign initial heads/elevation for nodes
  //! \param[in] node_head vector of pair of nodal index and initial nodal head
  void assign_node_head(const std::vector<std::pair<Index, double>>& node_head);
  void assign_node_elevation(const std::vector<std::pair<Index, double>>& node_head);

  //! Assign initial demand for nodes that have known discharge
  //! \param[in] node_discharge vector of pair of nodal index and initial
  //! discharge
  void assign_node_demand(
      const std::vector<std::pair<Index, double>>& node_discharge);

  //! Make MatrixAssembler a friend class of mesh
  friend class MatrixAssembler;

 private:
  //! mesh id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! nodal id and corresponding nodal pointer
  std::map<Index, std::shared_ptr<pipenetwork::Node>> nodes_;
  //! pipe id and corresponding pipe pointer
  std::map<Index, std::shared_ptr<pipenetwork::Pipe>> pipes_;
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_MESH_H_
