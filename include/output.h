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
#include <iomanip>

#include "mesh.h"
#ifndef PIPE_NETWORK_OUTPUT_H
#define PIPE_NETWORK_OUTPUT_H

namespace pipenetwork {
//! function converts input value to GPM unit for hydraulic simulation
//! \param[in] val the value need to be converted
//! \param[in] mode type of the input value (length, elevation, etc.)
double from_si(double val, const std::string&);

//! Pipe network input output class
//! \brief Base class for writing mesh information into .inp file
class Output {
 public:
  //! constructor
  //! \param[in] filename path of the .inp file
  Output(const std::shared_ptr<Mesh>& m, const std::string& filepath)
      : filepath_{filepath}, mesh_{m}, outFile_{filepath} {
    write_junctions();
    write_reservoirs();
    write_pipes();
  }

 private:
  //! Filepath
  std::string filepath_;
  //! Output object
  std::ofstream outFile_;
  //! Mesh pointer
  std::shared_ptr<Mesh> mesh_;
  //! Method for writing junctions
  void write_junctions();
  //! Method for writing reservoirs
  void write_reservoirs();
  //! Method for writing ouoes
  void write_pipes();
};
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_OUTPUT_H
