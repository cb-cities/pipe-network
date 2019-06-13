#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "settings.h"

#ifndef PIPE_NETWORK_INPUT_H
#define PIPE_NETWORK_INPUT_H

namespace pipenetwork {
//! Pipe network input output class
//! \brief Base class for parsing input data using .inp file
class Input {
  friend double to_si(double val, std::string mode);

 public:
  //! constructor
  //! \param[in] filename path of the .inp file
  explicit Input(const std::string & filename) : filename_(filename) {
    parse_sections();
    construct_node_info();
    construct_pipe_info();
  };

  //! Return node information

  std::vector<Index> junction_ids() const { return junction_ids_; }
  std::vector<double> junction_elevations() const {
    return junction_elevations_;
  }
  std::vector<double> junction_demands() const { return junction_demands_; }

  std::vector<Index> reservoir_ids() const { return reservoir_ids_; }
  std::vector<double> reservoir_heads() const { return reservoir_heads_; }

  //! Return pipe information
  std::vector<std::pair<Index, Index>> pipe_nodes_ids() const {
    return nodeids_;
  }
  std::vector<Index> pipe_ids() const { return pipe_ids_; }
  std::vector<double> pipe_diameters() const { return diameter_; }
  std::vector<double> pipe_length() const { return length_; }
  std::vector<double> pipe_roughness() const { return roughness_; }
  std::vector<Pipe_status> pipe_status() const { return pipe_status_; }

 private:
  //! Filepath
  std::string filename_;
  //! Set of section key words
  std::set<std::string> section_keys_{"[JUNCTIONS]", "[RESERVOIRS]", "[TANKS]",
                                      "[PIPES]", "[COORDINATES]"};
  //! number of attribute of the graph
  Index njunction_{0}, nresvoir_{0}, ntank_{0}, npipe_{0};

  //! Map of sections and its corresponding lines
  std::map<std::string, std::vector<std::string>> sections_;

  //! info for node constructions
  std::vector<std::pair<Index, Eigen::Vector3d>> node_coords_;
  std::vector<Index> junction_ids_;
  std::vector<double> junction_elevations_;
  std::vector<double> junction_demands_;

  std::vector<Index> reservoir_ids_;
  std::vector<double> reservoir_heads_;

  //! info for pipe constructions
  std::vector<Index> pipe_ids_;
  std::vector<std::pair<Index, Index>> nodeids_;
  std::vector<double> diameter_;
  std::vector<double> length_;
  std::vector<double> roughness_;
  std::vector<Pipe_status> pipe_status_;

  //! Parse information to each section from the input file
  void parse_sections();

  //! Construct node info that is used for pipeline-mesh
  void construct_node_info() {
    construct_node_elevation_head();
    construct_node_demand();
  }

  //  void construct_node_coord();
  void construct_node_elevation_head();
  void construct_node_demand();

  //! Construct pipe info that is used for pipeline-mesh
  void construct_pipe_info();

  //! Parse elevation, head or demand from a given node section
  //! \param[in] an opened file
  //!\retval a pair with vector of node ids their corresponding elevations,
  //!heads or demand
  std::pair<std::vector<Index>, std::vector<double>> parse_node_line(
      const std::string& section_name, const std::string& mode) const;
};

double to_si(double val, const std::string &);
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_INPUT_H
