#ifndef PIPE_NETWORK_IO_H
#define PIPE_NETWORK_IO_H

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
#include "io_utils.h"

namespace pipenetwork {

//! Pipe network input output class
//! \brief Base class for parsing input data using .inp file
class IO {
 public:
  //! Default constructor
  IO() = default;
  //! Read and create mesh from the .inp file
  //! \param[in] filename path of the .inp file
  void read_inp(const std::string& filename);

  //! Create synthetic network
  //! \param[in] n size of the network, number of nodes will be n^2
  void create_synthetic_net(Index n);

  //! save mesh information to a .inp file
  //! \param[in] mesh mesh ptr for save
  //! \param[in] output_path save path
  void save_mesh_inp(const std::shared_ptr<Mesh>& mesh,
                     const std::string& output_path);

  //! save hydrualic simulation result of the mesh
  //! \param[in] mesh mesh ptr for save
  //! \param[in] output_path save path
  void save_sim_result(const std::shared_ptr<Mesh>& mesh,
                       const std::string& output_path);

  //! Return node information
  std::vector<pipenetwork::JunctionProp> junction_properties() const {
    return junc_props_;
  }
  std::vector<pipenetwork::ReservoirProp> reservoir_properties() const {
    return res_props_;
  }
  //! Return pipe information
  std::vector<pipenetwork::PipeProp> pipe_properties() const {
    return pipe_props_;
  }
  //! Return pump information
  std::vector<PumpProp> pump_properties() const { return pump_props_; }

  //! Return valve information
  std::vector<ValveProp> valve_properties() const { return valve_props_; }

  //! Return curve information
  std::shared_ptr<Curves> curve_info() const { return curves_info_; }

 private:
  //! Filepath
  std::string filename_;
  //! the curves info ptr
  std::shared_ptr<Curves> curves_info_{std::make_shared<pipenetwork::Curves>()};
  //! Set of section key words
  std::set<std::string> section_keys_{
      "[JUNCTIONS]", "[RESERVOIRS]", "[TANKS]",  "[PIPES]", "[COORDINATES]",
      "[PUMPS]",     "[CURVES]",     "[VALVES]", "[LEAKS]"};

  //! Map of sections and its corresponding lines
  std::map<std::string, std::vector<std::string>> sections_;

  //! Map of leaking nid to leak diameter
  std::map<std::string, double> nid2ldia_;

  //! Vector of parsed junction properties
  std::vector<pipenetwork::JunctionProp> junc_props_;
  //! Vector of parsed reservoir properties
  std::vector<pipenetwork::ReservoirProp> res_props_;

  //! Vector of parsed pipe properties
  std::vector<PipeProp> pipe_props_;
  //! Vector of parsed pump properties
  std::vector<PumpProp> pump_props_;
  //! Vector of parsed valve properties
  std::vector<ValveProp> valve_props_;

  //! Parse information to each section from the input file
  void parse_sections();
  //! Parse leak information (leak diameter) from the input file
  void parse_leak_info();

  //! Parse elevation, head or demand from a given node section
  //! \param[in] an opened file
  //!\retval a pair with vector of node ids their corresponding elevations,
  //! heads or demand
  std::pair<std::vector<std::string>, std::vector<double>> parse_node_line(
      const std::string& section_name, const std::string& mode) const;

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
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_IO_H
