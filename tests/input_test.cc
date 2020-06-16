#include "catch.hpp"

#include "io.h"
#include "mesh.h"
// Check IO class
TEST_CASE("Input is checked", "[IO]") {
  double tolerance = 1e-6;

  SECTION("Check Parsed valve info") {

    // Create a INPUT class object
    auto IO = std::make_unique<pipenetwork::IO>();
    IO->read_inp("../test_files/test_net_valve.inp");

    // Mesh index
    std::string meshid = "IO_test_mesh";
    // Creat a mesh
    auto mesh = std::make_shared<pipenetwork::Mesh>(meshid);

    auto valve_props = IO->valve_properties();
    REQUIRE(valve_props.size() == 3);
    REQUIRE(valve_props[0].type == pipenetwork::ValveType::PRVALVE);
    REQUIRE(valve_props[0].name == "124");
    REQUIRE(valve_props[0].setting == 28.137549042234017);
    REQUIRE(valve_props[1].name == "125");
    REQUIRE(valve_props[1].type == pipenetwork::ValveType::TCVALVE);
    REQUIRE(valve_props[1].setting == 50);
  }

  SECTION("Check Parsed pump info") {
    // Create a INPUT class object
    auto IO = std::make_unique<pipenetwork::IO>();
    IO->read_inp("../test_files/test_net_pump.inp");
    // Mesh index
    std::string meshid = "IO_test_mesh";
    // Creat a mesh
    auto mesh = std::make_shared<pipenetwork::Mesh>(meshid);

    auto pump_props = IO->pump_properties();
    auto curve_info = IO->curve_info();
    REQUIRE(pump_props.size() == 2);
    REQUIRE(curve_info->pump_int_str(pump_props[0].curve_id) == "2");
    REQUIRE(pump_props[0].type == pipenetwork::PumpType::HEADPUMP);
    REQUIRE(pump_props[0].name == "123");
    REQUIRE(pump_props[1].name == "124");
    REQUIRE(pump_props[1].type == pipenetwork::PumpType::POWERPUMP);
    REQUIRE(pump_props[1].curve_id == -1);
    REQUIRE(pump_props[1].power == 37284.9936);
  }

  SECTION("Check Parsed node info") {
    auto IO = std::make_unique<pipenetwork::IO>();
    IO->read_inp("../test_files/test_net_valve.inp");
    // Mesh index
    std::string meshid = "IO_test_mesh";
    // Creat a mesh
    auto mesh = std::make_unique<pipenetwork::Mesh>(meshid);

    // junctions
    auto junction_props = IO->junction_properties();

    REQUIRE(junction_props.size() == 9);
    REQUIRE(junction_props[0].name == "10");
    REQUIRE(junction_props[0].elevation == Approx(216.408).epsilon(tolerance));
    REQUIRE(junction_props[1].demand ==
            Approx(0.0094635295).epsilon(tolerance));
    REQUIRE(junction_props[8].name == "32");
    REQUIRE(junction_props[0].leak_diameter == 0);

    // leak junction
    REQUIRE(junction_props[1].leak_diameter == 2 * 0.0254);

    // reservoirs
    auto reservoir_props = IO->reservoir_properties();

    REQUIRE(reservoir_props[0].head == Approx(243.84).epsilon(tolerance));
    REQUIRE(reservoir_props[0].name == "41");
  }
  SECTION("Check Parsed Pipe info") {
    auto IO = std::make_unique<pipenetwork::IO>();
    IO->read_inp("../test_files/test_net.inp");
    // check end nodes
    auto pipe_props = IO->pipe_properties();
    REQUIRE(pipe_props.size() == 8);
    REQUIRE(pipe_props[0].node1_name == "10");
    REQUIRE(pipe_props[0].node2_name == "11");
    // check roughness
    REQUIRE(pipe_props[0].roughness == 100);
    // check diameter
    REQUIRE(pipe_props[1].diameter == Approx(0.254).epsilon(tolerance));
    // check diameter
    REQUIRE(pipe_props[0].status == pipenetwork::LinkStatus::OPEN);
    // check length
    REQUIRE(pipe_props[1].length == Approx(1609.344).epsilon(tolerance));
    // check id
    REQUIRE(pipe_props[1].name == "102");
  }

  SECTION("Check Synthetic Net") {
    // Create a INPUT class object
    auto IO = std::make_shared<pipenetwork::IO>();
    IO->create_synthetic_net(10);
    auto junction_props = IO->junction_properties();
    REQUIRE(junction_props.size() == 100);
  }
}