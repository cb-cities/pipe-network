#include "catch.hpp"

#include "hydralic_sim.h"
// Check matrix_assembler class
TEST_CASE("HydraulicSimulation is checked", "[hydralic_sim]") {

  // Tolerance
  const double tolerance = 1.e-6;

  // Mesh index
  const unsigned meshid = 101;

  // Creat a mesh
  auto mesh = std::make_shared<pipenetwork::Mesh>(meshid);

  std::vector<Index> junction_ids{10, 11, 12};
  std::vector<double> elevations{216.408, 216.408, 213.36};
  std::vector<double> demands{0, 9.464e-03, 9.464e-03};
  std::vector<double> leak_diameters{0, 0.1, 0};

  mesh->create_junctions(junction_ids, elevations, demands, leak_diameters);

  std::vector<Index> res_ids{13};
  std::vector<double> heads{3.048e+02};

  mesh->create_reservoirs(res_ids, heads);

  std::vector<Index> pipe_ids{9, 10, 11, 12};
  std::vector<std::pair<Index, Index>> nodeids{
      std::make_pair(13, 10), std::make_pair(10, 11), std::make_pair(11, 12),
      std::make_pair(10, 12)};
  const std::vector<double> length{3209.5440000000003, 3209.5440000000003,
                                   1609.344, 1609.344};
  const std::vector<double> diameter{0.5588, 0.4572, 0.35559999999999997,
                                     0.254};
  const std::vector<double> roughness{100, 100, 100, 100};
  const std::vector<Pipe_status> status{OPEN, OPEN, OPEN, OPEN};

  mesh->create_pipes(pipe_ids, nodeids, length, diameter, roughness, status);

  double init_discharge = 1e-3;

  mesh->iterate_over_links(std::bind(&pipenetwork::Link::update_sim_discharge,
                                     std::placeholders::_1,
                                     init_discharge));  // initialze discharge
  SECTION("DD SIM TEST CASE 1: MESH INPUT") {
    bool pdd_mode = false;
    auto sim = std::make_shared<pipenetwork::Hydralic_sim>(mesh, pdd_mode);
    REQUIRE(sim->run_simulation());
    REQUIRE(sim->sim_residual_norm() < tolerance);

    // pressure test for stability
    REQUIRE(!sim->run_simulation(1e-30, 1000));
    REQUIRE(sim->sim_residual_norm() < tolerance);
  }

  SECTION("DD SIM TEST CASE 2: .INP FILE INPUT") {
    std::vector<double> leak_diameters{0, 0, 0, 0, 0, 0, 0, 0, 0};
    auto sim = std::make_shared<pipenetwork::Hydralic_sim>(
        "../benchmarks/Net1c.inp", leak_diameters);
    REQUIRE(sim->run_simulation());
    REQUIRE(sim->sim_residual_norm() < tolerance);
    // pressure test for stability
    REQUIRE(!sim->run_simulation(1e-30, 1000));
    REQUIRE(sim->sim_residual_norm() < tolerance);
  }
}
