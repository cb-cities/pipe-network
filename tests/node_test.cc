#include "catch.hpp"
#include <iostream>

#include "junction.h"
#include "reservoir.h"

// Check node class
TEST_CASE("Node is checked", "[Node]") {
  // Tolerance
  const double tolerance = 1.e-6;

  SECTION("Check junction") {
    pipenetwork::Junction_prop junction_1;
    junction_1.id = "123";
    junction_1.elevation = 10;
    junction_1.demand = 20;
    junction_1.leak_diameter = 5;

    auto junction_node = std::make_shared<pipenetwork::Junction>(junction_1);

    // check elevation
    REQUIRE(junction_node->id() == junction_1.id);

    // check elevation
    REQUIRE(junction_node->nodal_info()["elevation"] == junction_1.elevation);
    // check type
    REQUIRE(junction_node->nodal_info()["type"] == JUNCTION);
    // check demand
    REQUIRE(junction_node->nodal_info()["demand"] == junction_1.demand);
    // check leak hole area
    REQUIRE(junction_node->nodal_info()["leak_area"] ==
            Approx(19.634954084936).epsilon(tolerance));

    // set some simulation result
    junction_node->update_sim_head(99);
    REQUIRE(junction_node->sim_head() == Approx(99).epsilon(tolerance));
  }

  SECTION("Check reservoir") {
    // Index
    pipenetwork::Reservoir_prop res_1;
    res_1.id = "321";
    res_1.head = 10;

    auto reservoir_node = std::make_shared<pipenetwork::Reservoir>(res_1);

    // check elevation
    REQUIRE(reservoir_node->id() == res_1.id);
    // check elevation
    REQUIRE(reservoir_node->nodal_info()["head"] == res_1.head);
    // check type
    REQUIRE(reservoir_node->nodal_info()["type"] == RESERVOIR);
  }
}
