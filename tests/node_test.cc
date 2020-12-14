#include "catch.hpp"
#include <iostream>
#include <memory>

#include "junction.h"
#include "reservoir.h"

// Check node class
TEST_CASE("Node is checked", "[Node]") {
  using namespace pipenetwork;
  // Tolerance
  const double tolerance = 1.e-6;

  SECTION("Check junction") {
    // create junction
    JunctionProp junction1;
    junction1.name = "junc1";
    junction1.demand = 20;
    junction1.leak_diameter = 5;

    Index id = 1;

    auto junction_node = std::make_shared<Junction>(id, junction1);

    // check id
    REQUIRE(junction_node->id() == id);
    // check leak hole area
    REQUIRE(junction_node->leak_area() ==
            Approx(19.634954084936).epsilon(tolerance));

    // check properties
    auto prop = junction_node->property();
    // check name
    REQUIRE(prop.name == junction1.name);
    // check elevation
    REQUIRE(prop.elevation == junction1.elevation);
    // check demand
    REQUIRE(prop.demand == junction1.demand);
  }

  SECTION("Check reservoir") {
    // Index
    Index res_id = 1;
    // create reservoir
    ReservoirProp res1;

    res1.name = "res1";
    res1.head = 10;

    auto reservoir_node =
        std::make_shared<pipenetwork::Reservoir>(res_id, res1);

    // check id
    REQUIRE(reservoir_node->id() == res_id);
    // check properties
    auto prop = reservoir_node->property();
    // check name
    REQUIRE(prop.name == res1.name);
    // check elevation
    REQUIRE(prop.head == res1.head);

    // reset head
    reservoir_node->head(20);
    REQUIRE(reservoir_node->head() == 20);
  }
}
