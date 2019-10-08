#include "catch.hpp"

#include "hydralic_sim.h"
#include <chrono>

using namespace std::chrono;
// Check matrix_assembler class
TEST_CASE("HydraulicSimulation is checked", "[hydralic_sim]") {

  // Tolerance
  const double tolerance = 1.e-6;

  // Mesh index
  std::string meshid = "111";

  // Creat a mesh
  auto mesh = std::make_shared<pipenetwork::Mesh>(meshid);

  std::vector<std::string> junction_ids{"10", "11", "12"};
  std::vector<double> elevations{216.408, 216.408, 213.36};
  std::vector<double> demands{0, 9.464e-03, 9.464e-03};
  std::vector<double> leak_diameters{0, 0.1, 0};

  std::vector<pipenetwork::Junction_prop> junc_props;
  for (int i = 0; i < elevations.size(); ++i) {
    pipenetwork::Junction_prop junc_prop;
    junc_prop.id = junction_ids[i];
    junc_prop.elevation = elevations[i];
    junc_prop.demand = demands[i];
    junc_prop.leak_diameter = leak_diameters[i];

    junc_props.emplace_back(junc_prop);
  }

  mesh->create_junctions(junc_props);

  std::vector<std::string> res_ids{"13"};
  std::vector<double> heads{3.048e+02};
  std::vector<pipenetwork::Reservoir_prop> res_props;
  for (int i = 0; i < res_ids.size(); ++i) {
    pipenetwork::Reservoir_prop res_prop;
    res_prop.id = res_ids[i];
    res_prop.head = heads[i];

    res_props.emplace_back(res_prop);
  }
  mesh->create_reservoirs(res_props);

  mesh->create_reservoirs(res_props);

  std::vector<std::string> pipe_ids{"10", "11", "12", "13"};
  std::vector<std::pair<std::string, std::string>> nodeids{
      std::make_pair("13", "10"), std::make_pair("10", "11"),
      std::make_pair("11", "12"), std::make_pair("10", "12")};
  const std::vector<double> length{3209.5440000000003, 3209.5440000000003,
                                   1609.344, 1609.344};
  const std::vector<double> diameter{0.5588, 0.4572, 0.35559999999999997,
                                     0.254};
  const std::vector<double> roughness{100, 100, 100, 100};
  const std::vector<pipenetwork::Link_status> status{
      pipenetwork::OPEN, pipenetwork::OPEN, pipenetwork::OPEN,
      pipenetwork::OPEN};

  std::vector<pipenetwork::Pipe_prop> pipe_props;
  for (int i = 0; i < pipe_ids.size(); ++i) {
    pipenetwork::Pipe_prop pipe_prop;
    pipe_prop.id = pipe_ids[i];
    pipe_prop.length = length[i];
    pipe_prop.diameter = diameter[i];
    pipe_prop.roughness = roughness[i];
    pipe_prop.node1_id = nodeids[i].first;
    pipe_prop.node2_id = nodeids[i].second;
    pipe_prop.status = status[i];

    pipe_props.emplace_back(pipe_prop);
  }

  mesh->create_pipes(pipe_props);

  double init_discharge = 1e-3;

  mesh->iterate_over_links(std::bind(&pipenetwork::Link::update_sim_discharge,
                                     std::placeholders::_1,
                                     init_discharge));  // initialze discharge
  auto curves_info = std::make_shared<pipenetwork::Curves>();
  //  SECTION("DD SIM TEST CASE 1: MESH INPUT") {
  //    bool pdd_mode = false;
  //    bool debug = false;
  //    std::string solver_name = "mkl_pardiso";
  //    auto sim = std::make_shared<pipenetwork::Hydralic_sim>(
  //        mesh, curves_info, pdd_mode, solver_name, debug);
  //    REQUIRE(sim->run_simulation());
  //    REQUIRE(sim->sim_residual_norm() < tolerance);
  //    // pressure test for stability
  //    REQUIRE(!sim->run_simulation(1e-30, 1000));
  //    REQUIRE(sim->sim_residual_norm() < tolerance);
  //  }
  //
  //  SECTION("DD SIM TEST CASE 4: Large .INP FILE INPUT") {
  //    bool pdd_mode = false;
  //    bool debug = false;
  //    std::string mesh_name = "ky10";
  //    std::string solver_name = "mkl_pardiso";
  //    auto sim = std::make_shared<pipenetwork::Hydralic_sim>(
  //        "../benchmarks/ky10.inp", mesh_name, pdd_mode, solver_name, debug);
  //    auto start = high_resolution_clock::now();
  //
  //    sim->run_simulation(1e-8, 100);
  //    auto stop = high_resolution_clock::now();
  //    auto duration = duration_cast<seconds>(stop - start);
  //
  //    //    std::cout << duration.count() << std::endl;
  //    REQUIRE(sim->sim_residual_norm() < tolerance);
  //  }
  //
  //  SECTION(
  //      "DD SIM TEST CASE 7: Complete .INP FILE INPUT (with pump and valves)")
  //      {
  //    bool pdd_mode = false;
  //    bool debug = false;
  //    std::string mesh_name = "testnet";
  //    std::string solver_name = "mkl_pardiso";
  //    auto sim = std::make_shared<pipenetwork::Hydralic_sim>(
  //        "../benchmarks/test_net.inp", mesh_name, pdd_mode, solver_name,
  //        debug);
  //    auto start = high_resolution_clock::now();
  //
  //    //    REQUIRE(sim->run_simulation());
  //    sim->run_simulation(1e-8, 100);
  //    auto stop = high_resolution_clock::now();
  //    auto duration = duration_cast<milliseconds>(stop - start);
  //
  //    //    std::cout << duration.count() << std::endl;
  //    REQUIRE(sim->sim_residual_norm() < tolerance);
  //  }
  //
  //  SECTION("DD SIM TEST CASE 8: Broken .INP FILE INPUT (with pump and
  //  valves") {
  //    bool pdd_mode = false;
  //    bool debug = false;
  //    std::string mesh_name = "broken_net";
  //    std::string solver_name = "mkl_pardiso";
  //    auto sim = std::make_shared<pipenetwork::Hydralic_sim>(
  //        "../benchmarks/test_net_broken.inp", mesh_name, pdd_mode,
  //        solver_name, debug);
  //    auto start = high_resolution_clock::now();
  //
  //    //    REQUIRE(sim->run_simulation());
  //    sim->run_simulation(1e-8, 100);
  //    auto stop = high_resolution_clock::now();
  //    auto duration = duration_cast<milliseconds>(stop - start);
  //
  //    //    std::cout << duration.count() << std::endl;
  //    REQUIRE(sim->sim_residual_norm() < tolerance);
  //  }
  SECTION("DD SIM TEST CASE 9: Large Synthetic Network") {
    bool pdd_mode = false;
    bool debug = false;
    std::string solver_name = "cuda";
    auto start = high_resolution_clock::now();
    auto sim = std::make_shared<pipenetwork::Hydralic_sim>(20, pdd_mode,
                                                           solver_name, debug);
    //    REQUIRE(sim->run_simulation());
    sim->run_simulation(1e-8, 30);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    //      std::cout << duration.count() << std::endl;
    REQUIRE(sim->sim_residual_norm() < tolerance);
  }
  //  SECTION("DD SIM TEST CASE 9: Saved Synthetic Network") {
  //    bool pdd_mode = false;
  //    bool debug = false;
  //    std::string solver_name = "mkl_pardiso";
  //    auto start = high_resolution_clock::now();
  //    auto sim = std::make_shared<pipenetwork::Hydralic_sim>(
  //        "../benchmarks/synthetic_400.inp", "synthetic_400", pdd_mode,
  //        solver_name, debug);
  //    //    REQUIRE(sim->run_simulation());
  //    sim->run_simulation(1e-8, 100);
  //    auto stop = high_resolution_clock::now();
  //    auto duration = duration_cast<milliseconds>(stop - start);
  //
  //    std::cout << duration.count() << std::endl;
  //    REQUIRE(sim->sim_residual_norm() < tolerance);
  //  }

  //  SECTION("PDD SIM TEST CASE 1: Large real case .INP FILE INPUT") {
  //    bool pdd_mode = true;
  //    bool debug = false;
  //    std::string solver_name = "mkl_pardiso";
  //    auto sim = std::make_shared<pipenetwork::Hydralic_sim>(
  //        "../benchmarks/ky1.inp", "KY1_pdd", pdd_mode, solver_name, debug);
  //    REQUIRE(sim->run_simulation(1e-8, 30));
  //    REQUIRE(sim->sim_residual_norm() < tolerance);
  //  }
}
