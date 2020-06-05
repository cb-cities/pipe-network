#include "catch.hpp"
#include <iostream>
#include <string>

#include "junction.h"
#include "pipe.h"
#include "pump.h"
#include "reservoir.h"
#include "valve.h"

// Check node class
TEST_CASE("Link is checked", "[Link]") {
  // Tolerance
  const double tolerance = 1.e-6;

  pipenetwork::Index nid1 = 1;
  pipenetwork::Index nid2 = 2;

  pipenetwork::ReservoirProp res1;
  res1.name = "res1";
  res1.head = 10;

  pipenetwork::JunctionProp junc1;
  junc1.name = "123";
  junc1.elevation = 10;
  junc1.demand = 20;
  junc1.leak_diameter = 5;

  // Create a junction
  auto junction_node = pipenetwork::Junction(nid1, junc1);
  // Create a reservoir
  auto reservoir_node = pipenetwork::Reservoir(nid2, res1);

  SECTION("Check pipe") {
    pipenetwork::PipeProp pipe1;
    pipe1.name = "pipe1";
    pipe1.length = 10;
    pipe1.diameter = 8;
    pipe1.roughness = 0.5;
    pipenetwork::Index pipe_id = 1;

    auto pipe =
        pipenetwork::Pipe(pipe_id, junction_node, reservoir_node, pipe1);
    // check id
    REQUIRE(pipe.id() == pipe_id);

    // check two end points
    REQUIRE(pipe.nodes().first->id() == nid1);
    REQUIRE(pipe.nodes().second->id() == nid2);

    // check pipe property
    auto property = pipe.property();
    REQUIRE(property.name == pipe1.name);
    REQUIRE(property.length == pipe1.length);
  }

  SECTION("Check pump") {
    pipenetwork::PumpProp pump1;
    pump1.name = "pump1";
    pump1.status = pipenetwork::LinkStatus::CLOSED;
    pipenetwork::Index pump_id = 2;

    auto pump =
        pipenetwork::Pump(pump_id, junction_node, reservoir_node, pump1);
    // check id
    REQUIRE(pump.id() == pump_id);

    // check two end points
    REQUIRE(pump.nodes().first->id() == nid1);
    REQUIRE(pump.nodes().second->id() == nid2);

    // check pipe property
    auto property = pump.property();
    REQUIRE(property.name == pump1.name);
    REQUIRE(property.status == pump1.status);
  }

  SECTION("Check Valve") {
    pipenetwork::ValveProp valve1;
    valve1.name = "valve1";
    valve1.status = pipenetwork::LinkStatus::CLOSED;
    pipenetwork::Index pump_id = 3;

    auto valve =
        pipenetwork::Valve(pump_id, junction_node, reservoir_node, valve1);
    // check id
    REQUIRE(valve.id() == pump_id);

    // check two end points
    REQUIRE(valve.nodes().first->id() == nid1);
    REQUIRE(valve.nodes().second->id() == nid2);

    // check pipe property
    auto property = valve.property();
    REQUIRE(property.name == valve1.name);
    REQUIRE(property.status == valve1.status);
  }
}
