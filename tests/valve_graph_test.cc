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
  std::vector<std::string> pipe_names{"p1", "p2", "p3", "p4"};

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
  v1_prop.on_pipe = "p1";
  v1_prop.name = "v1";
  iso_valve_props.emplace_back(v1_prop);

  pipenetwork::isolation::ISOVProp v2_prop;
  v2_prop.on_node = "2";
  v2_prop.on_pipe = "p1";
  v2_prop.name = "v2";
  iso_valve_props.emplace_back(v2_prop);
  // assume broken

  pipenetwork::isolation::ISOVProp v3_prop;
  v3_prop.on_node = "2";
  v3_prop.on_pipe = "p3";
  v3_prop.name = "v3";
  iso_valve_props.emplace_back(v3_prop);

  pipenetwork::isolation::ISOVProp v4_prop;
  v4_prop.on_node = "3";
  v4_prop.on_pipe = "p3";
  v4_prop.name = "v4";
  iso_valve_props.emplace_back(v4_prop);

  pipenetwork::isolation::ISOVProp v5_prop;
  v5_prop.on_node = "4";
  v5_prop.on_pipe = "p4";
  v5_prop.name = "v5";
  iso_valve_props.emplace_back(v5_prop);

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
    //    std::cout << valve_loc_matrix << std::endl;
    REQUIRE(valve_loc_matrix.coeff(0, 0) == 1);
    REQUIRE(valve_loc_matrix.coeff(1, 0) == 1);
    REQUIRE(valve_loc_matrix.coeff(1, 1) == 0);
    REQUIRE(valve_loc_matrix.coeff(3, 1) == 0);
    REQUIRE(valve_loc_matrix.coeff(1, 2) == 1);
    REQUIRE(valve_loc_matrix.coeff(2, 2) == 1);
    REQUIRE(valve_loc_matrix.coeff(3, 2) == 0);
    REQUIRE(valve_loc_matrix.coeff(3, 3) == 1);
  }
  SECTION("CHECK isolation segment search algorithm") {
    auto& initial_segments = valve_graph.segments();

    auto nsegs = initial_segments.nsegs();
    REQUIRE(nsegs == 4);

    pipenetwork::Index pid = 0;
    auto seg0 = initial_segments.pid2seg(pid);

    REQUIRE(seg0.pids.size() == 1);
    REQUIRE(seg0.nids.size() == 0);
    REQUIRE(seg0.vids.size() == 2);

    auto seg3 = initial_segments.pid2seg(3);
    REQUIRE(seg3.pids.size() == 1);
    REQUIRE(seg3.nids.size() == 1);
    REQUIRE(seg3.vids.size() == 2);
  }
  SECTION("Check the Segment Valves Matrix") {
    auto& initial_segments = valve_graph.segments();
    auto seg_valve_mtx = initial_segments.seg_valve_mtx();
    std::cout << seg_valve_mtx << std::endl;
    REQUIRE(seg_valve_mtx.innerSize() == 4);
    REQUIRE(seg_valve_mtx.outerSize() == 5);

    REQUIRE(seg_valve_mtx.coeff(0, 0) == 1);
    REQUIRE(seg_valve_mtx.coeff(0, 1) == 1);
    REQUIRE(seg_valve_mtx.coeff(0, 2) == 0);
    REQUIRE(seg_valve_mtx.coeff(1, 0) == 0);
    REQUIRE(seg_valve_mtx.coeff(1, 1) == 1);
    REQUIRE(seg_valve_mtx.coeff(1, 2) == 1);
    REQUIRE(seg_valve_mtx.coeff(1, 4) == 1);
    REQUIRE(seg_valve_mtx.coeff(2, 0) == 0);
    REQUIRE(seg_valve_mtx.coeff(2, 2) == 1);
    REQUIRE(seg_valve_mtx.coeff(2, 3) == 1);
    REQUIRE(seg_valve_mtx.coeff(3, 1) == 0);
    REQUIRE(seg_valve_mtx.coeff(3, 3) == 1);
    REQUIRE(seg_valve_mtx.coeff(3, 4) == 1);

    SECTION("Check the shrink matrix tool") {
      std::vector<pipenetwork::Index> rows_to_remove{0, 2};
      std::vector<pipenetwork::Index> cols_to_remove{1, 3};
      auto shrinked_mtx = pipenetwork::isolation::IsoSegHelper::shrink_mtx(
          seg_valve_mtx, rows_to_remove, cols_to_remove);
      REQUIRE(shrinked_mtx.innerSize() == 2);
      REQUIRE(shrinked_mtx.outerSize() == 3);
      REQUIRE(shrinked_mtx.coeff(0, 0) == 0);
      REQUIRE(shrinked_mtx.coeff(0, 1) == 1);
      REQUIRE(shrinked_mtx.coeff(0, 2) == 1);
      REQUIRE(shrinked_mtx.coeff(1, 0) == 0);
      REQUIRE(shrinked_mtx.coeff(1, 1) == 0);
      REQUIRE(shrinked_mtx.coeff(1, 2) == 1);
    }
  }
  SECTION("Check the Segment Valves ADJ Matrix") {
    auto& initial_segments = valve_graph.segments();
    auto seg_valve_adj_mtx = initial_segments.seg_valve_adj_mtx();
    //      std::cout << seg_valve_adj_mtx << std::endl;
    REQUIRE(seg_valve_adj_mtx.coeff(0, 0) == 0);
    REQUIRE(seg_valve_adj_mtx.coeff(0, 1) == 1);
    REQUIRE(seg_valve_adj_mtx.coeff(1, 0) == 1);
    REQUIRE(seg_valve_adj_mtx.coeff(1, 1) == 0);
  }
  SECTION("Check components finding algorithm") {
    auto& segments = valve_graph.segments();
    std::set<pipenetwork::Index> iso_sid{1};
    auto components = segments.get_segment_components(iso_sid);
    REQUIRE(components.size() == 3);
    REQUIRE(components[0].size() == 1);
    REQUIRE(components[1].size() == 1);
    REQUIRE(components[2].size() == 2);
  }
  SECTION("Check merge segments by remove valves") {
    auto& segments = valve_graph.segments();
    segments.merge_segments(1);
    auto seg_valve_mtx = segments.seg_valve_mtx();
    std::cout << seg_valve_mtx << std::endl;
    REQUIRE(seg_valve_mtx.innerSize() == 3);
    REQUIRE(seg_valve_mtx.coeff(0, 0) == 1);
    REQUIRE(seg_valve_mtx.coeff(0, 1) == 1);
    REQUIRE(seg_valve_mtx.coeff(0, 2) == 1);
    REQUIRE(seg_valve_mtx.coeff(0, 4) == 1);
    REQUIRE(seg_valve_mtx.coeff(1, 0) == 0);
    REQUIRE(seg_valve_mtx.coeff(1, 1) == 0);
    REQUIRE(seg_valve_mtx.coeff(1, 2) == 1);
    REQUIRE(seg_valve_mtx.coeff(1, 4) == 0);
    REQUIRE(seg_valve_mtx.coeff(2, 0) == 0);
    REQUIRE(seg_valve_mtx.coeff(2, 1) == 0);
    REQUIRE(seg_valve_mtx.coeff(2, 2) == 0);
    REQUIRE(seg_valve_mtx.coeff(2, 4) == 1);

    auto segment = segments.pid2seg(1);
    REQUIRE(segment.sid == 0);
    REQUIRE(segment.pids.size() == 2);

    SECTION("Check components finding algorithm after merging") {
      std::set<pipenetwork::Index> iso_sid{2};
      auto components = segments.get_segment_components(iso_sid);
      REQUIRE(components.size() == 2);
      REQUIRE(components[0].size() == 2);
      REQUIRE(components[1].size() == 1);
    }
    SECTION("Continous failure case") {
      segments.merge_segments(4);
      auto seg_valve_mtx = segments.seg_valve_mtx();
      std::cout << seg_valve_mtx << std::endl;
      REQUIRE(seg_valve_mtx.innerSize() == 2);

      auto seg_valve_adj_mtx = segments.seg_valve_adj_mtx();
      std::cout << seg_valve_adj_mtx << std::endl;

      auto components = segments.get_segment_components({1});
      REQUIRE(components.size() == 2);
      REQUIRE(components[0].size() == 1);
      REQUIRE(components[1].size() == 1);
    }
  }

  SECTION("Check wrappers in the valve graph class") {
    valve_graph.initialize_iso_segs();  // reinitialize
    REQUIRE(valve_graph.nsegs() == 4);

    auto seg0 = valve_graph.pname2segment("p1");
    REQUIRE(seg0.pids.size() == 1);
    REQUIRE(seg0.nids.size() == 0);
    REQUIRE(seg0.vids.size() == 2);

    auto seg_components_init = valve_graph.segment_components({"p1", "p2"});
    REQUIRE(seg_components_init.size() == 3);

    valve_graph.fail_valves({"v2", "v5"});
    REQUIRE(valve_graph.nsegs() == 2);

    auto seg_components_after = valve_graph.segment_components({"p1", "p2"});
    REQUIRE(seg_components_after.size() == 2);

    auto seg3 = valve_graph.pname2segment("p3");
    REQUIRE(seg3.pids.size() == 1);
    REQUIRE(seg3.nids.size() == 0);
    REQUIRE(seg3.vids.size() == 2);

    auto seg_union = valve_graph.pname2segment("p4");
    REQUIRE(seg_union.pids.size() == 3);
    REQUIRE(seg_union.nids.size() == 3);
    REQUIRE(seg_union.vids.size() == 3);
  }
}