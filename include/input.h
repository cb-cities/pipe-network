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
  Input(std::string filename) : filename_(filename) {
    parse_sections();
    construct_node_info();
    construct_pipe_info();
  };

  //! Return node information
  std::vector<std::pair<Index, Eigen::Vector3d>> get_node_coord() const {
    return node_coords_;
  }
  std::vector<std::pair<Index, double>> get_node_elevation() const {
    return node_elevations_;
  }
  std::vector<std::pair<Index, double>> get_node_demand() const {
    return node_demands_;
  }

  //! Return pipe information
  std::vector<std::pair<Index, Index>> get_pipe_end_ids() const {
    return nodeids_;
  }
  std::vector<double> get_pipe_diameters() const { return diameter_; }
  std::vector<double> get_pipe_length() const { return length_; }
  std::vector<double> get_pipe_roughness() const { return roughness_; }
  std::vector<bool> get_pipe_status() const { return pipe_status_; }

  //! Return node id map
  //! \retval a map with node id in file and its corresponding id in computation
  //! mesh
  std::map<Index, Index> get_node_map() { return node_id_map_; }

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

  //! Map of node id between file and mesh
  std::map<Index, Index> node_id_map_;
  //! Map of pipe id between file and mesh
  std::map<Index, Index> pipe_id_map_;

  //! info for node constructions
  std::vector<std::pair<Index, Eigen::Vector3d>> node_coords_;
  std::vector<std::pair<Index, double>> node_elevations_;
  std::vector<std::pair<Index, double>> node_demands_;

  //! info for pipe constructions
  std::vector<std::pair<Index, Index>> nodeids_;
  std::vector<double> diameter_;
  std::vector<double> length_;
  std::vector<double> roughness_;
  std::vector<bool> pipe_status_;

  //! Parse information to each section from the input file
  void parse_sections();

  //! Construct node info that is used for pipeline-mesh
  void construct_node_info() {
    construct_node_coord();
    construct_node_elevation();
    construct_node_demand();
  }
  void construct_node_coord();
  void construct_node_elevation();
  void construct_node_demand();

  //! Construct pipe info that is used for pipeline-mesh
  void construct_pipe_info();

  //! Parse elevation or demand from a given node section
  //! \param[in] an opened file
  //!\retval a vector with node ids in mesh and their corresponding elevations
  std::vector<std::pair<Index, double>> parse_node_line(
      std::string section_name, std::string mode) const;
};

double to_si(double val, std::string mode);
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_INPUT_H
