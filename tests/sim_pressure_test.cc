#include "catch.hpp"
#include "hydralic_sim.h"
#include <chrono>

using namespace std::chrono;
// Check matrix_assembler class
TEST_CASE("HydraulicSimulation pressure test", "[hydralic_sim]") {

  // Tolerance
  const double tolerance = 1.e-6;

  SECTION("DD SIM TEST CASE 4: Large .INP FILE INPUT") {
    bool pdd_mode = false;
    bool debug = false;
    std::string mesh_name = "ky7";
    std::string solver_name = "mkl_pardiso";
    auto sim = std::make_shared<pipenetwork::Hydralic_sim>(
        "../benchmarks/ky7.inp", mesh_name, pdd_mode, solver_name, debug);
    auto start = high_resolution_clock::now();

    sim->run_simulation(1e-8, 100);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);

    //    std::cout << duration.count() << std::endl;
    REQUIRE(sim->sim_residual_norm() < tolerance);
  }
  SECTION("DD SIM TEST CASE 8: Broken .INP FILE INPUT (with pump and valves") {
    bool pdd_mode = false;
    bool debug = false;
    std::string mesh_name = "broken_net";
    std::string solver_name = "mkl_pardiso";
    auto sim = std::make_shared<pipenetwork::Hydralic_sim>(
        "../test_files/test_net_broken.inp", mesh_name, pdd_mode, solver_name,
        debug);
    auto start = high_resolution_clock::now();

    //    REQUIRE(sim->run_simulation());
    sim->run_simulation(1e-8, 100);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    //      std::cout << duration.count() << std::endl;
    REQUIRE(sim->sim_residual_norm() < tolerance);
  }
  SECTION("DD SIM TEST CASE 9: Large Synthetic Network") {
    bool pdd_mode = false;
    bool debug = false;
    std::string solver_name = "mkl_pardiso";
    auto start = high_resolution_clock::now();
    auto sim = std::make_shared<pipenetwork::Hydralic_sim>(100, pdd_mode,
                                                           solver_name, debug);
    //    REQUIRE(sim->run_simulation());
    sim->run_simulation(1e-8, 30);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    //    std::cout << duration.count() << std::endl;
    REQUIRE(sim->sim_residual_norm() < tolerance);
  }

  SECTION("PDD SIM TEST CASE 1: Large real case .INP FILE INPUT") {
    bool pdd_mode = true;
    bool debug = false;
    std::string solver_name = "mkl_pardiso";
    auto sim = std::make_shared<pipenetwork::Hydralic_sim>(
        "../benchmarks/ky7.inp", "KY7_pdd", pdd_mode, solver_name, debug);
    REQUIRE(sim->run_simulation(1e-8, 30));
    REQUIRE(sim->sim_residual_norm() < tolerance);
  }
}
