#include "catch.hpp"

#include "hydralic_sim.h"
#include <chrono>

using namespace std::chrono;
// Check matrix_assembler class
TEST_CASE("HydraulicSimulation is checked", "[hydralic_sim]") {

  // Tolerance
  const double tolerance = 1.e-8;

  // Mesh index
  std::string meshid = "Matrix test mesh";

  // Creat a mesh
  auto mesh = std::make_shared<pipenetwork::Mesh>(meshid);

  // Create a curves info object
  auto curves_info = std::make_shared<pipenetwork::Curves>();

  std::vector<std::string> junction_ids{"10", "11", "12"};
  std::vector<double> elevations{216.408, 216.408, 213.36};
  std::vector<double> demands{0, 9.464e-03, 9.464e-03};
  std::vector<double> leak_diameters{0, 0.1, 0};

  std::vector<pipenetwork::JunctionProp> junc_props;
  for (int i = 0; i < elevations.size(); ++i) {
    pipenetwork::JunctionProp junc_prop;
    junc_prop.name = junction_ids[i];
    junc_prop.elevation = elevations[i];
    junc_prop.demand = demands[i];
    junc_prop.leak_diameter = leak_diameters[i];

    junc_props.emplace_back(junc_prop);
  }

  std::vector<std::string> res_ids{"13"};
  std::vector<double> heads{3.048e+02};

  std::vector<pipenetwork::ReservoirProp> res_props;
  for (int i = 0; i < res_ids.size(); ++i) {
    pipenetwork::ReservoirProp res_prop;
    res_prop.name = res_ids[i];
    res_prop.head = heads[i];

    res_props.emplace_back(res_prop);
  }

  std::vector<std::string> pipe_ids{"10", "11"};
  std::vector<std::pair<std::string, std::string>> nodeids{
      std::make_pair("13", "10"), std::make_pair("10", "11")};
  const std::vector<double> length{3209.5440000000003, 3209.5440000000003,
                                   1609.344, 1609.344};
  const std::vector<double> diameter{0.5588, 0.4572, 0.35559999999999997,
                                     0.254};
  const std::vector<double> roughness{100, 100, 100, 100};

  std::vector<pipenetwork::PipeProp> pipe_props;
  for (int i = 0; i < pipe_ids.size(); ++i) {
    pipenetwork::PipeProp pipe_prop;
    pipe_prop.name = pipe_ids[i];
    pipe_prop.length = length[i];
    pipe_prop.diameter = diameter[i];
    pipe_prop.roughness = roughness[i];
    pipe_prop.node1_name = nodeids[i].first;
    pipe_prop.node2_name = nodeids[i].second;

    pipe_props.emplace_back(pipe_prop);
  }

  std::vector<std::string> pump_ids{"12"};
  std::vector<std::pair<std::string, std::string>> pump_nodeids{
      std::make_pair("10", "12")};

  std::vector<pipenetwork::PumpProp> pump_props;

  for (int i = 0; i < pump_ids.size(); ++i) {
    pipenetwork::PumpProp pump_prop;
    pump_prop.name = pump_ids[i];
    pump_prop.type = pipenetwork::PumpType::POWERPUMP;
    pump_prop.node1_name = pump_nodeids[i].first;
    pump_prop.node2_name = pump_nodeids[i].second;
    pump_props.emplace_back(pump_prop);
  }

  std::vector<std::string> valve_ids{"13"};
  std::vector<std::pair<std::string, std::string>> valve_nodeids{
      std::make_pair("11", "12")};

  std::vector<pipenetwork::ValveProp> valve_props;

  for (int i = 0; i < valve_ids.size(); ++i) {
    pipenetwork::ValveProp valve_prop;
    valve_prop.name = valve_ids[i];
    valve_prop.type = pipenetwork::ValveType::PRVALVE;
    valve_prop.status = pipenetwork::LinkStatus ::ACTIVE;
    valve_prop.setting = 10;
    valve_prop.node1_name = valve_nodeids[i].first;
    valve_prop.node2_name = valve_nodeids[i].second;

    valve_props.emplace_back(valve_prop);
  }

  mesh->create_nodes(junc_props, res_props);
  mesh->create_links(pipe_props, pump_props, valve_props);
  mesh->create_mesh_graph();

  SECTION("DD SIM TEST CASE 1: MESH INPUT ") {
    bool pdd_mode = false;
    bool debug = false;
    auto sim = std::make_shared<pipenetwork::Hydralic_sim>(mesh, curves_info,
                                                           pdd_mode, debug);
    REQUIRE(sim->run_simulation());
    REQUIRE(sim->sim_residual_norm() < tolerance);
    // pressure test for stability
    auto start = high_resolution_clock::now();
    REQUIRE(!sim->run_simulation(1e-30, 100));
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    //    std::cout << duration.count() << std::endl;
    REQUIRE(sim->sim_residual_norm() < tolerance);
    sim->update_mesh();
  }

  SECTION("DD SIM TEST CASE 2: .INP FILE INPUT") {
    bool pdd_mode = false;
    bool debug = false;

    auto IO = std::make_shared<pipenetwork::IO>();
    IO->read_inp("../test_files/test_net.inp");
    std::string mesh_name = "test_mesh_inp";

    auto mesh_inp = std::make_shared<pipenetwork::Mesh>(mesh_name);
    mesh_inp->create_nodes(IO->junction_properties(),
                           IO->reservoir_properties());
    mesh_inp->create_links(IO->pipe_properties(), IO->pump_properties(),
                           IO->valve_properties());
    mesh_inp->create_mesh_graph();
    mesh_inp->print_summary();

    auto curves_info_io = IO->curve_info();
    auto sim = std::make_shared<pipenetwork::Hydralic_sim>(
        mesh_inp, curves_info_io, pdd_mode, debug);
    sim->run_simulation(1e-8, 10);
    REQUIRE(sim->sim_residual_norm() < tolerance);
  }

  SECTION("PDD SIM TEST CASE: .INP FILE INPUT") {
    bool pdd_mode = true;
    bool debug = false;

    auto IO = std::make_shared<pipenetwork::IO>();
    IO->read_inp("../test_files/test_net_pdd.inp");
    std::string mesh_name = "pdd_test";

    auto mesh_inp = std::make_shared<pipenetwork::Mesh>(mesh_name);
    mesh_inp->create_nodes(IO->junction_properties(),
                           IO->reservoir_properties());
    mesh_inp->create_links(IO->pipe_properties(), IO->pump_properties(),
                           IO->valve_properties());
    mesh_inp->create_mesh_graph();
    mesh_inp->print_summary();

    auto curves_info_io = IO->curve_info();
    auto sim = std::make_shared<pipenetwork::Hydralic_sim>(
        mesh_inp, curves_info_io, pdd_mode, debug);
    sim->run_simulation(1e-8, 20);
    REQUIRE(sim->sim_residual_norm() < tolerance);
  }
}
