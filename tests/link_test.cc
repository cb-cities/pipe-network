#include "catch.hpp"
#include <iostream>
#include <string>

#include "junction.h"
#include "pipe.h"
#include "reservoir.h"

// Check node class
TEST_CASE("Link is checked", "[Link]") {
  // Tolerance
  const double tolerance = 1.e-6;

  pipenetwork::Reservoir_prop res_1;
  res_1.id = "321";
  res_1.head = 10;

  pipenetwork::Junction_prop junction_1;
  junction_1.id = "123";
  junction_1.elevation = 10;
  junction_1.demand = 20;
  junction_1.leak_diameter = 5;

  auto junction_node = std::make_shared<pipenetwork::Junction>(junction_1);
  // Create a reservoir
  auto reservoir_node = std::make_shared<pipenetwork::Reservoir>(res_1);

  SECTION("Check pipe") {
    pipenetwork::Pipe_prop pipe1;
    pipe1.id = "111";
    pipe1.length = 10;
    pipe1.diameter = 8;
    pipe1.roughness = 0.5;
    pipe1.status = pipenetwork::CLOSED;
    pipe1.node1 = reservoir_node;
    pipe1.node2 = junction_node;
    auto pipe = std::make_shared<pipenetwork::Pipe>(pipe1);
    // check id
    REQUIRE(pipe->id() == pipe1.id);

    // check two end points
    REQUIRE(pipe->nodes().first->id() == res_1.id);
    REQUIRE(pipe->nodes().second->id() == junction_1.id);

    // check pipe type
    REQUIRE(pipe->link_info()["type"] == pipenetwork::PIPE);
    // check length
    REQUIRE(pipe->link_info()["length"] ==
            Approx(pipe1.length).epsilon(tolerance));
    // check diameter
    REQUIRE(pipe->link_info()["diameter"] ==
            Approx(pipe1.diameter).epsilon(tolerance));
    // check roughness
    REQUIRE(pipe->link_info()["roughness"] ==
            Approx(pipe1.roughness).epsilon(tolerance));
    // check status
    REQUIRE(pipe->link_status() == pipenetwork::CLOSED);
  }
}
