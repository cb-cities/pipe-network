#include "catch.hpp"

#include "input.h"
#include "mesh.h"
// Check IO class
TEST_CASE("Input is checked", "[IO]") {
  double tolerance = 1e-6;
  // Create a INPUT class object
  auto IO = std::make_shared<pipenetwork::Input>("../benchmarks/test_net.inp");
  // Mesh index
  const unsigned meshid = 101;
  // Creat a mesh
  auto mesh = std::make_unique<pipenetwork::Mesh>(meshid);
  SECTION("Check Parsed pump info") {
    auto valve_props = IO->valve_properties();
    REQUIRE(valve_props.size() == 2);
    REQUIRE(valve_props[0].valve_type == pipenetwork::PRVALVE);
    REQUIRE(valve_props[0].id == "124");
    REQUIRE(valve_props[0].setting == 28.137549042234017);
    REQUIRE(valve_props[1].id == "125");
    REQUIRE(valve_props[1].valve_type == pipenetwork::TCVALVE);
    REQUIRE(valve_props[1].setting == 50);

    mesh->create_mesh_from_inp(IO);
    mesh->print_summary();
  }

  //  SECTION("Check Parsed pump info") {
  //    auto pump_props = IO->pump_properties();
  //    auto curve_info = IO->curve_info();
  //    REQUIRE(pump_props.size() == 2);
  //    REQUIRE(curve_info->pump_int_str(pump_props[0].curve_name) == "2");
  //    REQUIRE(pump_props[0].pump_type == pipenetwork::HEADPUMP);
  //    REQUIRE(pump_props[0].id == "123");
  //    REQUIRE(pump_props[1].id == "124");
  //    REQUIRE(pump_props[1].pump_type == pipenetwork::POWERPUMP);
  //    REQUIRE(pump_props[1].curve_name == -1);
  //    REQUIRE(pump_props[1].power == 37284.9936);
  //    auto junction_props = IO->junction_properties();
  //    auto reservoir_props = IO->reservoir_properties();
  //
  //    mesh->create_junctions (junction_props);
  //    mesh->create_reservoirs (reservoir_props);
  //
  //    mesh->create_pumps (pump_props);
  ////    mesh->print_summary ();
  //  }
  //  SECTION("Check Parsed node info") {
  //
  //    // junctions
  //    auto junction_props = IO->junction_properties();
  //
  //    REQUIRE(junction_props.size() == 9);
  //    REQUIRE(junction_props[0].id == "10");
  //    REQUIRE(junction_props[0].elevation ==
  //    Approx(216.408).epsilon(tolerance)); REQUIRE(junction_props[1].demand ==
  //            Approx(0.0094635295).epsilon(tolerance));
  //    REQUIRE(junction_props[8].id == "32");
  //
  //    // reservoirs
  //    auto reservoir_props = IO->reservoir_properties();
  //
  //    REQUIRE(reservoir_props[0].head == Approx(304.8).epsilon(tolerance));
  //    REQUIRE(reservoir_props[0].id == "99");
  //  }
  //  SECTION("Check Parsed Pipe info") {
  //    // check end nodes
  //    auto pipe_props = IO->pipe_properties();
  //    REQUIRE(pipe_props.size() == 13);
  //    REQUIRE(pipe_props[0].node1_id == "10");
  //    REQUIRE(pipe_props[0].node2_id == "11");
  //    // check roughness
  //    REQUIRE(pipe_props[0].roughness == 100);
  //    // check diameter
  //    REQUIRE(pipe_props[1].diameter ==
  //            Approx(0.35559999999999997).epsilon(tolerance));
  //    // check diameter
  //    REQUIRE(pipe_props[0].status == pipenetwork::OPEN);
  //    // check length
  //    REQUIRE(pipe_props[1].length == Approx(1609.344).epsilon(tolerance));
  //    // check id
  //    REQUIRE(pipe_props[12].id == "110");
  //  }
}