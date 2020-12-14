#ifndef PIPE_NETWORK_IO_UTILS_H
#define PIPE_NETWORK_IO_UTILS_H

#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

#include "mesh.h"

namespace pipenetwork {
namespace IO_utils {
//! function to convert input value to standard unit for hydraulic simulation
//! \param[in] val the value need to be converted
//! \param[in] mode type of the input value (length, elevation, etc.)
double to_si(double val, const std::string&);

//! function converts input value to GPM unit for hydraulic simulation
//! \param[in] val the value need to be converted
//! \param[in] mode type of the input value (length, elevation, etc.)
double from_si(double val, const std::string&);

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
//! function check the existance of a path
//! \param[in] pathname folder path
bool path_exist(const std::string& pathname);

//! function to create a new folder if it does not exist
//! \param[in] pathname folder path
void create_new_folder(const std::string& pathname);

//! Pipe network input output class
//! \brief Base class for writing mesh information into .inp file
class Output {
 public:
  //! constructor
  //! \param[in] m mesh ptr
  explicit Output(const std::shared_ptr<Mesh>& m) : mesh_{m} {}

  //! save mesh information to a .inp file
  //! \param[in] output_path path for save
  void save_mesh_inp(const std::string& output_path);

  //! save hydraulic simulation result
  //! \param[in] output_path path for save
  void save_sim_result(const std::string& output_path);

 private:
  //! Filepath
  std::string filepath_;
  //! Mesh pointer
  std::shared_ptr<Mesh> mesh_;

  //! Method for writing junctions
  void write_junctions_inp(std::ofstream& outfile);

  //! Method for writing reservoirs
  void write_reservoirs_inp(std::ofstream& outfile);

  //! Method for writing ouoes
  void write_pipes_inp(std::ofstream& outfile);
};
}  // namespace IO_utils
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_IO_UTILS_H
