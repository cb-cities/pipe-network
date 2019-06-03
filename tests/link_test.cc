#include "catch.hpp"
#include <iostream>

#include "pipe.h"
#include "junction.h"
#include "reservoir.h"

// Check node class
TEST_CASE("Link is checked", "[Link]") {
    // Tolerance
    const double tolerance = 1.e-6;
    // Index
    unsigned junction_id{333},reservoir_id{222};
    // Create a junction
    double elevation = 10;
    double demand = 20;
    double leak_diameter = 5;
    auto junction_node = std::make_shared<pipenetwork::Junction>(
            junction_id, elevation, demand, leak_diameter);
    //Create a reservoir
    double head = 999;
    auto reservoir_node = std::make_shared<pipenetwork::Reservoir>(reservoir_id, head);

    SECTION("Check pipe") {
        Index pipe_id = 111;
        double length{10},diameter{8},roughness{0.5};
        auto pipe = std::make_shared<pipenetwork::Pipe>(pipe_id, reservoir_node,
                junction_node,length,diameter,roughness,closed);
        // check id
        REQUIRE(pipe->id() == pipe_id);

        // check two end points
        REQUIRE(pipe->nodes ().first->id () == reservoir_id);
        REQUIRE(pipe->nodes ().second->id () == junction_id);

        // check pipe type
        REQUIRE(pipe->link_info ()["type"] == pipe_type);
        // check length
        REQUIRE(pipe->link_info()["length"] == Approx(length).epsilon(tolerance));
        // check diameter
        REQUIRE(pipe->link_info()["diameter"] == Approx(diameter).epsilon(tolerance));
        // check roughness
        REQUIRE(pipe->link_info()["roughness"] == Approx(roughness).epsilon(tolerance));
        // check status
        REQUIRE(pipe->link_info()["status"] == closed);
    }

}

