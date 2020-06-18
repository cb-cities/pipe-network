#include <chrono>
#include <exception>
#include <string>
#include <tclap/CmdLine.h>

#include "hydralic_sim.h"

int main(int argc, char** argv) {
  using namespace std::chrono;
  try {
    TCLAP::CmdLine cmd(
        "HydroSim: fast hydraulic simulation for potable water distribution "
        "system",
        ' ', "1.0");
    TCLAP::ValueArg<std::string> fileArg(
        "f", "file", ".inp file for WDN system", false,
        "../test_files/test_net.inp", "filename", cmd);
    TCLAP::ValueArg<std::string> saveArg("t", "to",
                                         "Folder to save the results", false,
                                         "../results/", "save path", cmd);
    TCLAP::ValueArg<std::string> nameArg("n", "name", "Name for the WDN", false,
                                         "test", "meshname", cmd);
    TCLAP::SwitchArg pddSwitch(
        "p", "pdd", "Run pressure demand driven (PDD) simulation", cmd, false);
    TCLAP::SwitchArg debugSwitch(
        "d", "debug", "Debug mode: output intermediate results", cmd, false);

    // Parse the argv array.
    cmd.parse(argc, argv);
    std::string filepath = fileArg.getValue();
    std::string mesh_name = nameArg.getValue();
    bool pdd_mode = pddSwitch.getValue();
    bool debug = debugSwitch.getValue();
    std::string save_path = saveArg.getValue();
    // IO
    auto IO = std::make_shared<pipenetwork::IO>();
    IO->read_inp(filepath);
    auto start = high_resolution_clock::now();
    // Mesh
    auto mesh = std::make_shared<pipenetwork::Mesh>(mesh_name);
    mesh->create_nodes(IO->junction_properties(), IO->reservoir_properties());
    mesh->create_links(IO->pipe_properties(), IO->pump_properties(),
                       IO->valve_properties());
    mesh->create_mesh_graph();
    mesh->print_summary();
    // hydraulic simulation
    auto curves_info_io = IO->curve_info();
    auto sim = std::make_shared<pipenetwork::Hydralic_sim>(mesh, curves_info_io,
                                                           pdd_mode, debug);
    if (sim->run_simulation(1e-8, 100)) {
      std::cout << "Simulation Completed!!!" << std::endl;
      sim->update_mesh();
      IO->save_sim_result(mesh, save_path);
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stop - start);
      std::cout << "hydraulic computation time (milliseconds)"
                << duration.count() << std::endl;
    } else {
      std::cerr << "Simulation does not converge!  " << std::endl;
    }
  } catch (std::exception& e) {
    std::cerr << "Simulation failed: " << e.what() << std::endl;
    std::abort();
  }

  return 0;
}
