#include "catch.hpp"

#include "input.h"
// Check IO class
TEST_CASE("Input is checked", "[IO]") {
  double tolerance = 1e-6;
  // Create a IO class object
  auto IO = std::make_unique<pipenetwork::Input>("../benchmarks/Net1c.inp");
  SECTION("Check Parsed node info") {

    // junctions
    auto junction_props = IO->junction_properties();

    REQUIRE(junction_props.size() == 9);
    REQUIRE(junction_props[0].id == "10");
    REQUIRE(junction_props[0].elevation == Approx(216.408).epsilon(tolerance));
    REQUIRE(junction_props[1].demand ==
            Approx(0.0094635295).epsilon(tolerance));
    REQUIRE(junction_props[8].id == "32");

    // reservoirs
    auto reservoir_props = IO->reservoir_properties();

    REQUIRE(reservoir_props[0].head == Approx(304.8).epsilon(tolerance));
    REQUIRE(reservoir_props[0].id == "99");
  }
  SECTION("Check Parsed Pipe info") {
    // check end nodes
    auto pipe_props = IO->pipe_properties();
    REQUIRE(pipe_props.size() == 13);
    REQUIRE(pipe_props[0].node1_id == "10");
    REQUIRE(pipe_props[0].node2_id == "11");
    // check roughness
    REQUIRE(pipe_props[0].roughness == 100);
    // check diameter
    REQUIRE(pipe_props[1].diameter ==
            Approx(0.35559999999999997).epsilon(tolerance));
    // check diameter
    REQUIRE(pipe_props[0].status == pipenetwork::OPEN);
    // check length
    REQUIRE(pipe_props[1].length == Approx(1609.344).epsilon(tolerance));
    // check id
    REQUIRE(pipe_props[12].id == "110");
  }
}