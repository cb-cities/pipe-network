#include <exception>
#include <string>
#include <tclap/CmdLine.h>

#include "hydralic_sim.h"

int main(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd(
        "HydroSim: fast hydraulic simulation for potable water distribution "
        "system",
        ' ', "1.0");
    TCLAP::ValueArg<std::string> fileArg(
        "f", "file", ".inp file for WDN system", false, "../benchmarks/ky7.inp",
        "filename", cmd);
    TCLAP::ValueArg<std::string> saveArg("t", "to",
                                         "Folder to save the results", false,
                                         "../results/", "save path", cmd);
    TCLAP::ValueArg<std::string> nameArg("n", "name", "Name for the WDN", false,
                                         "test", "meshname", cmd);
    TCLAP::ValueArg<std::string> solverArg("s", "solver", "Solver to use",
                                           false, "mkl_pardiso", "solvername",
                                           cmd);
    TCLAP::SwitchArg pddSwitch(
        "p", "pdd", "Run pressure demand driven (PDD) simulation", cmd, false);
    TCLAP::SwitchArg debugSwitch(
        "d", "debug", "Debug mode: output intermediate results", cmd, false);

    // Parse the argv array.
    cmd.parse(argc, argv);

    std::string filepath = fileArg.getValue();
    std::string mesh_name = nameArg.getValue();
    std::string solver_name = solverArg.getValue();
    bool pdd_mode = pddSwitch.getValue();
    bool debug = debugSwitch.getValue();
    std::string save_path = saveArg.getValue();

    auto sim = std::make_shared<pipenetwork::Hydralic_sim>(
        filepath, mesh_name, pdd_mode, solver_name, debug);
    if (sim->run_simulation(1e-8, 30, save_path)) {
      std::cout << "Simulation Completed!!!" << std::endl;
    } else {
      std::cerr << "Simulation does not converge!  " << std::endl;
    }

  } catch (std::exception& e) {
    std::cerr << "Simulation failed: " << e.what() << std::endl;
  }

  return 0;
}
