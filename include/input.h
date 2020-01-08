#ifndef PIPE_NETWORK_INPUT_H
#define PIPE_NETWORK_INPUT_H

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
#include "valve.h"

namespace pipenetwork {
//! function to convert input value to standard unit for hydraulic simulation
//! \param[in] val the value need to be converted
//! \param[in] mode type of the input value (length, elevation, etc.)
double to_si(double val, const std::string&);

//! function to generate uniformly distributed random number
//! \param[in] l min value boundary for the random number domain
//! \param[in] r max value boundary for the random number domain
double rand_number(double l, double h);

//! function check the existance of a file
//! \param[in] name path of the file
inline bool file_exists(const std::string& name) {
  std::ifstream f(name.c_str());
  return f.good();
}

//! Pipe network input output class
//! \brief Base class for parsing input data using .inp file
class Input {
 public:
  //! constructor to create synthetic mesh
  //! \param[in] n number of the mesh dimension (n*n)
  explicit Input(int n) {
    auto junction_nodes = construct_synthesis_junctions(n);
    construct_synthesis_pipes(junction_nodes);
    create_sources(n);
  };

  //! constructor to create mesh from .inp file
  //! \param[in] filename path of the .inp file
  explicit Input(const std::string& filename) : filename_(filename) {
    if (!file_exists(filename)) {
      throw std::runtime_error(
          "Input file does not exist, please check the path!");
    }
    parse_sections();
    construct_node_info();
    construct_pipe_info();
    construct_curve_info();
    construct_pump_info();
    construct_valve_info();
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

  //! Return valve information
  std::vector<Valve_prop> valve_properties() const { return valve_props_; }

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
                                      "[CURVES]",    "[VALVES]"};

  //! Map of sections and its corresponding lines
  std::map<std::string, std::vector<std::string>> sections_;

  //! Vector of parsed junction properties
  std::vector<pipenetwork::Junction_prop> junc_props_;
  //! Vector of parsed reservoir properties
  std::vector<pipenetwork::Reservoir_prop> res_props_;

  //! Vector of parsed pipe properties
  std::vector<Pipe_prop> pipe_props_;
  //! Vector of parsed pump properties
  std::vector<Pump_prop> pump_props_;
  //! Vector of parsed valve properties
  std::vector<Valve_prop> valve_props_;

  //! Parse information to each section from the input file
  void parse_sections();

  //! Construct node info that is used for pipeline-mesh
  void construct_node_info();
  //! Construct pipe info that is used for pipeline-mesh
  void construct_pipe_info();
  //! Construct curve info that is used for pipeline-mesh and matrix assembler
  void construct_curve_info();
  //! Construct pump info that is used for pipeline-mesh
  void construct_pump_info();
  //! Construct valve info that is used for pipeline-mesh
  void construct_valve_info();

  //! Construct synthesis junctions
  //! \param[in] n number of the mesh dimension (n*n)
  std::vector<std::vector<std::string>> construct_synthesis_junctions(int n);

  //! Construct synthesis pipes
  //! \param[in] junction_names names for junctions
  void construct_synthesis_pipes(
      const std::vector<std::vector<std::string>>& junction_names);
  //! Construct pipes by linking junctions vertically
  //! \param[in] junction_names names for junctions
  //! \param[in] col_num current column number
  //! \param[in] rand random connection or not
  void create_vertical_pipes(const std::vector<std::string>& junction_names,
                             int col_num, bool rand = false);
  //! Construct pipes by linking junctions horizontally
  //! \param[in] l_junc names for left sided junctions
  //! \param[in] r_junc names for right sided junctions
  //! \param[in] col_num current column number
  //! \param[in] rand random connection or not
  void create_horizontal_pipes(const std::vector<std::string>& l_junc,
                               const std::vector<std::string>& r_junc,
                               int col_num, bool rand = true);
  //! Create sources for synthetic network
  //! \param[in] n number of sources for the synthetic network
  void create_sources(int n);

  //! Parse elevation, head or demand from a given node section
  //! \param[in] an opened file
  //!\retval a pair with vector of node ids their corresponding elevations,
  //! heads or demand
  std::pair<std::vector<std::string>, std::vector<double>> parse_node_line(
      const std::string& section_name, const std::string& mode) const;
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_INPUT_H
