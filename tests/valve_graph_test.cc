#include "catch.hpp"

#include "io.h"
#include "valve_graph.h"

// Check mesh class
TEST_CASE("Valve Graph is checked", "[Valve Graph]") {
  // Tolerance
  const double tolerance = 1.e-12;

  // create junctions
  std::vector<std::string> junction_names{"1", "2", "3", "4", "5"};
  std::vector<double> elevations{5, 4, 3, 2, 1};
  std::vector<double> demands{1, 2, 3, 4, 5};
  std::vector<double> leak_diameters{0, 0, 0, 0.2, 0.1};

  std::vector<pipenetwork::JunctionProp> junc_props;
  for (int i = 0; i < elevations.size(); ++i) {
    pipenetwork::JunctionProp junc_prop;
    junc_prop.name = junction_names[i];
    junc_prop.elevation = elevations[i];
    junc_prop.demand = demands[i];
    junc_prop.leak_diameter = leak_diameters[i];
    junc_props.emplace_back(junc_prop);
  }

  // Reservoirs
  std::vector<std::string> res_names{"6", "7"};
  std::vector<double> heads{99, 100};

  std::vector<pipenetwork::ReservoirProp> res_props;
  for (int i = 0; i < res_names.size(); ++i) {
    pipenetwork::ReservoirProp res_prop;
    res_prop.name = res_names[i];
    res_prop.head = heads[i];
    res_props.emplace_back(res_prop);
  }

  std::vector<std::string> pipe_names{"1", "2", "3"};
  // 1<->2<->6; 4<->7

  std::vector<std::pair<std::string, std::string>> node_names{
      std::make_pair("1", "2"), std::make_pair("2", "6"),
      std::make_pair("4", "7")};
  const std::vector<double> length{100, 200, 300};
  const std::vector<double> diameter{3, 4, 5};
  const std::vector<double> roughness{0.2, .6, .9};

  std::vector<pipenetwork::PipeProp> pipe_props;
  for (int i = 0; i < pipe_names.size(); ++i) {
    pipenetwork::PipeProp pipe_prop;
    pipe_prop.name = pipe_names[i];
    pipe_prop.length = length[i];
    pipe_prop.diameter = diameter[i];
    pipe_prop.roughness = roughness[i];
    pipe_prop.node1_name = node_names[i].first;
    pipe_prop.node2_name = node_names[i].second;
    pipe_props.emplace_back(pipe_prop);
  }

  std::vector<pipenetwork::PumpProp> pump_props;
  std::vector<pipenetwork::ValveProp> valve_props;

  for (int i = 0; i < 3; ++i) {
    pipenetwork::PumpProp pump_prop;
    pump_prop.name = "pump-" + std::to_string(i);
    pump_props.emplace_back(pump_prop);
  }

  std::vector<pipenetwork::ISOVProp> iso_valve_props;
  pipenetwork::ISOVProp v1_prop;
  v1_prop.on_node = "2";
  v1_prop.on_pipe = "2";
  iso_valve_props.emplace_back(v1_prop);

  pipenetwork::ISOVProp v2_prop;
  v2_prop.on_node = "4";
  v2_prop.on_pipe = "3";
  iso_valve_props.emplace_back(v2_prop);

  auto mesh_nodes =
      std::make_shared<pipenetwork::MeshNodes>(junc_props, res_props);
  auto mesh_links = std::make_shared<pipenetwork::MeshLinks>(
      pipe_props, pump_props, valve_props, *mesh_nodes);
  auto valve_graph =
      pipenetwork::ValveGraph(mesh_nodes, mesh_links, iso_valve_props);

  SECTION("CHECK The NODE PIPE MATRIX") {
    auto node_pipe_matrix = valve_graph.node_pipe_mtx();
    REQUIRE(node_pipe_matrix.coeff(0, 0) == 1);
    REQUIRE(node_pipe_matrix.coeff(1, 0) == 1);
    REQUIRE(node_pipe_matrix.coeff(1, 1) == 1);
    REQUIRE(node_pipe_matrix.coeff(5, 1) == 1);
    REQUIRE(node_pipe_matrix.coeff(3, 2) == 1);
    REQUIRE(node_pipe_matrix.coeff(6, 2) == 1);
    REQUIRE(node_pipe_matrix.coeff(5, 2) == 0);
  }
  SECTION("CHECK The Valve Location Matrix") {
    auto valve_loc_matrix = valve_graph.valve_loc_mtx();
    REQUIRE(valve_loc_matrix.coeff(1, 1) == 1);
    REQUIRE(valve_loc_matrix.coeff(3, 2) == 1);
    REQUIRE(valve_loc_matrix.coeff(5, 1) == 0);
    REQUIRE(valve_loc_matrix.coeff(6, 2) == 0);
    REQUIRE(valve_loc_matrix.coeff(5, 2) == 0);
  }
  SECTION("CHECK isolation segment search algorithm") {
      pipenetwork::Index pid = 0;
      auto seg0 = valve_graph.get_iso_seg (pid);

      REQUIRE(seg0.pids.size() == 1);
//      REQUIRE(seg0.pids[0] == 0);

      REQUIRE(seg0.nids.size() == 2);
//      REQUIRE(seg0.nids[0] == 1);
//      REQUIRE(seg0.nids[1] == 0);
  }
}