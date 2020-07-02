#include "catch.hpp"

#include "io.h"
#include "valve_graph.h"
#include <Eigen/Eigenvalues>

// Check mesh class
TEST_CASE("Valve Graph is checked", "[Valve Graph]") {
  // Tolerance
  const double tolerance = 1.e-12;

  // create junctions
  std::vector<std::string> junction_names{"1", "2", "3", "4"};
  std::vector<double> elevations{5, 4, 3, 2};
  std::vector<double> demands{1, 2, 3, 4};
  std::vector<double> leak_diameters{0, 0, 0, 0.2};

  std::vector<pipenetwork::JunctionProp> junc_props;
  for (int i = 0; i < elevations.size(); ++i) {
    pipenetwork::JunctionProp junc_prop;
    junc_prop.name = junction_names[i];
    junc_prop.elevation = elevations[i];
    junc_prop.demand = demands[i];
    junc_prop.leak_diameter = leak_diameters[i];
    junc_props.emplace_back(junc_prop);
  }

  std::vector<pipenetwork::ReservoirProp> res_props;
  std::vector<std::string> pipe_names{"1", "2", "3", "4"};

  std::vector<std::pair<std::string, std::string>> node_names{
      std::make_pair("1", "2"), std::make_pair("2", "4"),
      std::make_pair("2", "3"), std::make_pair("3", "4")};
  const std::vector<double> length{100, 200, 300, 400};
  const std::vector<double> diameter{3, 4, 5, 6};
  const std::vector<double> roughness{0.2, .6, .9, 2};

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

  std::vector<pipenetwork::isolation::ISOVProp> iso_valve_props;
  pipenetwork::isolation::ISOVProp v1_prop;
  v1_prop.on_node = "1";
  v1_prop.on_pipe = "1";
  v1_prop.name = "v1";
  iso_valve_props.emplace_back(v1_prop);

  pipenetwork::isolation::ISOVProp v2_prop;
  v2_prop.on_node = "2";
  v2_prop.on_pipe = "1";
  v2_prop.name = "v2";
  // assume broken

  pipenetwork::isolation::ISOVProp v3_prop;
  v3_prop.on_node = "2";
  v3_prop.on_pipe = "3";
  v3_prop.name = "v3";
  iso_valve_props.emplace_back(v3_prop);

  pipenetwork::isolation::ISOVProp v4_prop;
  v4_prop.on_node = "3";
  v4_prop.on_pipe = "3";
  v4_prop.name = "v4";
  iso_valve_props.emplace_back(v4_prop);

  pipenetwork::isolation::ISOVProp v5_prop;
  v5_prop.on_node = "4";
  v5_prop.on_pipe = "2";
  v5_prop.name = "v5";
  iso_valve_props.emplace_back(v5_prop);

  pipenetwork::isolation::ISOVProp v6_prop;
  v6_prop.on_node = "4";
  v6_prop.on_pipe = "4";
  v6_prop.name = "v6";
  iso_valve_props.emplace_back(v6_prop);

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
    REQUIRE(node_pipe_matrix.coeff(3, 1) == 1);
    REQUIRE(node_pipe_matrix.coeff(1, 2) == 1);
    REQUIRE(node_pipe_matrix.coeff(2, 2) == 1);
    REQUIRE(node_pipe_matrix.coeff(5, 2) == 0);
  }
  SECTION("CHECK The Valve Location Matrix") {
    auto valve_loc_matrix = valve_graph.valve_loc_mtx();
    REQUIRE(valve_loc_matrix.coeff(0, 0) == 1);
    REQUIRE(valve_loc_matrix.coeff(1, 0) == 0);  // assume broken
    REQUIRE(valve_loc_matrix.coeff(1, 1) == 0);
    REQUIRE(valve_loc_matrix.coeff(3, 1) == 1);
    REQUIRE(valve_loc_matrix.coeff(1, 2) == 1);
    REQUIRE(valve_loc_matrix.coeff(2, 2) == 1);
    REQUIRE(valve_loc_matrix.coeff(5, 2) == 0);
    //      std::cout<<valve_loc_matrix<<std::endl;
  }
  SECTION("CHECK isolation segment search algorithm") {

    auto nsegs = valve_graph.nsegs();
    REQUIRE(nsegs == 3);

    pipenetwork::Index pid = 0;
    auto seg0 = valve_graph.get_iso_seg(pid);

    REQUIRE(seg0.pids.size() == 2);
    REQUIRE(seg0.nids.size() == 1);
    REQUIRE(seg0.vids.size() == 3);
  }
  SECTION("Check the Segment Valves Matrix") {
    auto seg_valve_mtx = valve_graph.seg_valve_mtx();

    REQUIRE(seg_valve_mtx.coeff(0, 0) == 1);
    REQUIRE(seg_valve_mtx.coeff(0, 1) == 1);
    REQUIRE(seg_valve_mtx.coeff(0, 2) == 0);
    REQUIRE(seg_valve_mtx.coeff(1, 0) == 0);
    REQUIRE(seg_valve_mtx.coeff(1, 2) == 1);

    std::cout << seg_valve_mtx << std::endl;

    auto new_seg_valve = valve_graph.merge_segments(1);
    std::cout << new_seg_valve << std::endl;
    //      seg_valve_mtx.col(2) += 3 * seg_valve_mtx.col(0);
    //      std::cout << seg_valve_mtx << std::endl;
    //    auto dense =   Eigen::MatrixXd(seg_valve_mtx);
    //    pipenetwork::isolation::IsoSegHelper::removeRow(dense, 0);
    //    std::cout
    //        << seg_valve_mtx << std::endl;
  }

  SECTION("Check the Segment Valves ADJ Matrix") {
    auto seg_valve_adj_mtx = valve_graph.seg_valve_adj_mtx();
    std::cout << seg_valve_adj_mtx << std::endl;
    REQUIRE(seg_valve_adj_mtx.coeff(0, 0) == 0);
    REQUIRE(seg_valve_adj_mtx.coeff(0, 1) == 1);
    REQUIRE(seg_valve_adj_mtx.coeff(1, 0) == 1);
    REQUIRE(seg_valve_adj_mtx.coeff(1, 1) == 0);
    valve_graph.find_segment_components(1);
  }
}