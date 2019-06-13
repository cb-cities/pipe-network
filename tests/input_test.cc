#include "catch.hpp"

#include "input.h"
// Check IO class
TEST_CASE("Input is checked", "[IO]") {
    double tolerance = 1e-6;
    // Create a IO class object
    auto IO = std::make_unique<pipenetwork::Input>("../benchmarks/Net1c.inp");
    SECTION ("Check Parsed node info"){

        // junctions
        auto junction_ids = IO->junction_ids ();
        auto junction_elevations = IO->junction_elevations ();
        auto junction_demands = IO->junction_demands ();

        REQUIRE (junction_ids.size() == 9);
        REQUIRE (junction_ids[0] == 10);
        REQUIRE (junction_elevations[0] ==  Approx(216.408).epsilon (tolerance));
        REQUIRE (junction_demands[1] == Approx(0.0094635295).epsilon (tolerance));
        REQUIRE (junction_ids[8] == 32);

        // reservoirs
        auto reservoir_id = IO->reservoir_ids ();
        auto reservoir_head = IO->reservoir_heads ();

        REQUIRE (reservoir_head[0] == Approx(304.8).epsilon (tolerance));
        REQUIRE (reservoir_id[0] == Approx(99).epsilon (tolerance));

    }
    SECTION ("Check Parsed Pipe info"){
        // check end nodes
        auto pipe_end_nodes = IO->pipe_nodes_ids ();
        REQUIRE (pipe_end_nodes.size() == 13);
        REQUIRE (pipe_end_nodes[0].first == 10);
        REQUIRE (pipe_end_nodes[0].second== 11);
        // check roughness
        auto rough = IO->pipe_roughness ();
        REQUIRE (rough[0] == 100);
        // check diameter
        auto dia = IO->pipe_diameters ();
        REQUIRE (dia[1] == Approx(0.35559999999999997).epsilon (tolerance));
        // check diameter
        auto status = IO->pipe_status ();
        REQUIRE (status[0] == OPEN);
        // check length
        auto length = IO->pipe_length ();
        REQUIRE (length[1] == Approx(1609.344).epsilon (tolerance));
        // check id
        REQUIRE (IO->pipe_ids ()[12] ==110);
    }





}