#ifndef PIPE_NETWORK_MESH_H
#define PIPE_NETWORK_MESH_H

#include <array>
#include <cmath>
#include <exception>
#include <iostream>
#include <map>
#include <memory>
#include <tuple>
#include <vector>

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
  bool create_junctions(const std::vector<Index>& ids,
                        const std::vector<double>& elevations,
                        const std::vector<double>& demands,
                        const std::vector<double>& leak_diameters);

  //! Create Reservoir pointers
  //! \param[in] ids reservoir id
  //! \param[in] heads base head for reservoirs
  bool create_reservoirs(const std::vector<Index>& ids,
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

  //! Print summary for the mesh
  void print_summary();

 private:
  //! mesh id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! nodal id and corresponding nodal pointer
  std::map<Index, std::shared_ptr<pipenetwork::Node>> nodes_;
  //! pipe id and corresponding pipe pointer
  std::map<Index, std::shared_ptr<pipenetwork::Link>> pipes_;
  //! nodal id and corresponding nodal pointer
  std::map<Index, std::shared_ptr<pipenetwork::Node>> connected_nodes_;
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_MESH_H
