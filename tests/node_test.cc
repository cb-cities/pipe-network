#include "catch.hpp"
#include <iostream>

#include "junction.h"
#include "reservoir.h"

// Check node class
TEST_CASE("Node is checked", "[Node]") {
  // Tolerance
  const double tolerance = 1.e-6;

  SECTION("Check junction") {
    // Index
    unsigned id = 123;
    double elevation = 10;
    double demand = 20;
    double leak_diameter = 5;
    auto junction_node = std::make_shared<pipenetwork::Junction>(
        id, elevation, demand, leak_diameter);

    // check elevation
    REQUIRE(junction_node->id() == id);

    // check elevation
    REQUIRE(junction_node->nodal_info()["elevation"] == elevation);
    // check type
    REQUIRE(junction_node->nodal_info()["type"] == junction_type);
    // check demand
    REQUIRE(junction_node->nodal_info()["demand"] == demand);
    // check leak hole area
    REQUIRE(junction_node->nodal_info()["leak_area"] ==
            Approx(19.634954084936).epsilon(tolerance));

    // set some simulation result
    junction_node->sim_head(99);
    REQUIRE(junction_node->sim_head() == Approx(99).epsilon(tolerance));
  }

  SECTION("Check reservoir") {
    // Index
    unsigned id = 321;
    double head = 10;
    auto reservoir_node = std::make_shared<pipenetwork::Reservoir>(id, head);

    // check elevation
    REQUIRE(reservoir_node->id() == id);
    // check elevation
    REQUIRE(reservoir_node->nodal_info()["head"] == head);
    // check type
    REQUIRE(reservoir_node->nodal_info()["type"] == reservoir_type);
  }
}
