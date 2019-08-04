#include <Eigen/Dense>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "curves.h"
#include "junction.h"
#include "pipe.h"
#include "pump.h"
#include "reservoir.h"
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
  explicit Input(const std::string& filename) : filename_(filename) {
    parse_sections();
    construct_node_info();
    construct_pipe_info();
    construct_curve_info();
    construct_pump_info();
  };

  //! Return node information
  std::vector<pipenetwork::Junction_prop> junction_properties() const {
    return junc_props_;
  }
  std::vector<pipenetwork::Reservoir_prop> reservoir_properties() const {
    return res_props_;
  }
  //! Return pipe information
  std::vector<pipenetwork::Pipe_prop> pipe_properties() const {
    return pipe_props_;
  }
  //! Return pump information
  std::vector<Pump_prop> pump_properties() const { return pump_props_; }

  //! Return curve information
  std::shared_ptr<Curves> curve_info() const { return curves_info_; }

 private:
  //! Filepath
  std::string filename_;
  //! the curves info ptr
  std::shared_ptr<Curves> curves_info_{std::make_shared<pipenetwork::Curves>()};
  //! Set of section key words
  std::set<std::string> section_keys_{"[JUNCTIONS]", "[RESERVOIRS]",  "[TANKS]",
                                      "[PIPES]",     "[COORDINATES]", "[PUMPS]",
                                      "[CURVES]"};
  //! number of attribute of the graph
  Index njunction_{0}, nresvoir_{0}, ntank_{0}, npipe_{0};

  //! Map of sections and its corresponding lines
  std::map<std::string, std::vector<std::string>> sections_;

  //! info for node constructions
  std::vector<std::pair<std::string, Eigen::Vector3d>> node_coords_;
  std::vector<std::string> junction_ids_;
  std::vector<double> junction_elevations_;
  std::vector<double> junction_demands_;
  std::vector<pipenetwork::Junction_prop> junc_props_;

  std::vector<std::string> reservoir_ids_;
  std::vector<double> reservoir_heads_;
  std::vector<pipenetwork::Reservoir_prop> res_props_;

  //! info for pipe constructions
  std::vector<std::string> pipe_ids_;
  std::vector<std::pair<std::string, std::string>> nodeids_;
  std::vector<double> diameter_;
  std::vector<double> length_;
  std::vector<double> roughness_;
  std::vector<Pipe_status> pipe_status_;
  std::vector<Pipe_prop> pipe_props_;
  std::vector<Pump_prop> pump_props_;

  //! Parse information to each section from the input file
  void parse_sections();

  //! Construct node info that is used for pipeline-mesh
  void construct_node_info();

  //  void construct_node_coord();
  void construct_node_elevation_head();
  void construct_node_demand();

  //! Construct pipe info that is used for pipeline-mesh
  void construct_pipe_info();
  void construct_curve_info();
  void construct_pump_info();

  //! Parse elevation, head or demand from a given node section
  //! \param[in] an opened file
  //!\retval a pair with vector of node ids their corresponding elevations,
  //! heads or demand
  std::pair<std::vector<std::string>, std::vector<double>> parse_node_line(
      const std::string& section_name, const std::string& mode) const;
};

double to_si(double val, const std::string&);
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_INPUT_H
