#ifndef PIPE_NETWORK_IO_H_
#define PIPE_NETWORK_IO_H_

#include <exception>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "csv.h"
#include <Eigen/Sparse>

#include "settings.h"

namespace pipenetwork {

//! Pipe network input output class
//! \brief Base class for input and output data
class IO {
 public:
  //! Read input CSV file and store the data
  //! \param[in] node_filename the name of input CSV file with nodal information
  //! \param[in] pipe_filename the name of input CSV file with pipe information
  bool read_network(const std::string& node_filename,
                    const std::string& pipe_filename);

  //! Return input nodal coordinates
  //! \retval nodal_coords_ input nodal coordinates
  std::vector<Eigen::Vector3d> nodal_coordinates() const {
    return nodal_coords_;
  }

  //! Return input start and end nodes of pipes
  //! \retval node_pairs input start and end nodes of pipes
  std::vector<std::pair<Index, Index>> node_pairs() const {
    return node_pairs_;
  }

  //! Return input pipe diameters
  //! \retval diameters_ input pipe diameters
  std::vector<double> diameters() const { return diameters_; }

  //! Return input pipe roughnesses
  //! \retval roughness_ input pipe diameters
  std::vector<double> roughness() const { return roughness_; }

  //! Return input pipe status
  //! \retval pipe_status_ input pipe status
  std::vector<bool> pipe_status() const { return pipe_status_; }

  //! Return input pipe discharge and corresponding index
  //! \retval initial_pipe_discharge_ input pipe discharge and corresponding
  //! index
  std::vector<std::pair<Index, double>> initial_pipe_discharge() const {
    return initial_pipe_discharge_;
  }

  //! Return input nodal head and corresponding index
  //! \retval initial_nodal_head_ input nodal head and corresponding index
  std::vector<std::pair<Index, double>> initial_nodal_head() const {
    return initial_nodal_head_;
  }

  //! Return input nodal discharge and corresponding index
  //! \retval initial_nodal_discharge_ input nodal discharge and corresponding
  //! index
  std::vector<std::pair<Index, double>> initial_nodal_discharge() const {
    return initial_nodal_discharge_;
  }

 private:
  //! Vector of nodal coordinates
  std::vector<Eigen::Vector3d> nodal_coords_;
  //! Vecotr of start and end nodes of pipes
  std::vector<std::pair<Index, Index>> node_pairs_;
  //! vector of pipe diameters
  std::vector<double> diameters_;
  //! vector of pipe roughnesses
  std::vector<double> roughness_;
  //! vector of pipe status
  std::vector<bool> pipe_status_;
  //! vector of pipe discharge and corresponding index
  std::vector<std::pair<Index, double>> initial_pipe_discharge_;
  //! vector of nodal head and corresponding index
  std::vector<std::pair<Index, double>> initial_nodal_head_;
  //! vector of nodal discharge and corresponding index
  std::vector<std::pair<Index, double>> initial_nodal_discharge_;
};
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_IO_H_
