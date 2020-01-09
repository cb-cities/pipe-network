#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "hydralic_sim.h"

bool to_bool(std::string str) {
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  std::istringstream is(str);
  bool b = (str == "true") || (str == "t");
  return b;
}

int main(int argc, char** argv) {
  std::string filepath = argv[1];
  std::string mesh_name = argv[2];
  std::string solver_name = argv[3];
  bool pdd_mode = to_bool(argv[4]);
  bool debug = to_bool(argv[5]);
  std::string save_path = argv[6];

  auto sim = std::make_shared<pipenetwork::Hydralic_sim>(
      filepath, mesh_name, pdd_mode, solver_name, debug);
  sim->run_simulation(1e-8, 30,save_path);
}
