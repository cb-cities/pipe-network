#include "catch.hpp"

#include "input.h"


// Check IO class
TEST_CASE("Input is checked", "[IO]") {
    double tolerance = 1e-6;
    // Create a IO class object
    auto IO = std::make_unique<pipenetwork::Input>("../benchmarks/Net1c.inp");
    SECTION ("Check Parsed node info"){
        auto node_info = IO->get_node_coord ();
        Index node_id = node_info[0].first;
        Eigen::Vector3d coord = node_info[0].second;


        REQUIRE (node_info.size() == 11);//check size
        REQUIRE (node_id == 0);// check ID
        // check coords
        REQUIRE (coord[0]== Approx(20).epsilon(tolerance));
        REQUIRE (coord[1]== Approx(70).epsilon(tolerance));
        REQUIRE (coord[2]== Approx(0).epsilon(tolerance));

//        for (auto const& x : IO->get_node_map ()){
//            std::cout<< x.first<<" "<<x.second<<std::endl;
//        }
        // node elevations
        auto node_elevation = IO->get_node_elevation ();
        //check size
        REQUIRE (node_elevation.size() == 11);
        //check a junction elevation
        REQUIRE (node_elevation[0].second == Approx(216.408).epsilon (tolerance));
//        //check a reservoir elevation
//        REQUIRE (node_elevation[92].second == Approx(220).epsilon (tolerance));
//        //check a tank elevation
//        REQUIRE (node_elevation[96].second == Approx(129).epsilon (tolerance));

        // node demand
        auto node_demand = IO->get_node_demand ();
        //check size
        REQUIRE (node_demand.size() == 11);
        //check a junction elevation
        REQUIRE (node_demand[1].second == Approx(0.0094635295).epsilon (tolerance));
//        //check a reservoir elevation
//        REQUIRE (node_demand[92].second == Approx(-1).epsilon (tolerance));
//        //check a tank elevation
//        REQUIRE (node_demand[96].second == Approx(-1).epsilon (tolerance));




    }
    SECTION ("Check Parsed Pipe info"){
        // check end nodes
        auto pipe_end_nodes = IO->get_pipe_end_ids ();
        REQUIRE (pipe_end_nodes.size() == 13);
        REQUIRE (pipe_end_nodes[0].first == 0);
        REQUIRE (pipe_end_nodes[0].second== 1);
        // check roughness
        auto rough = IO->get_pipe_roughness ();
        REQUIRE (rough[0] == 100);
        // check diameter
        auto dia = IO->get_pipe_diameters ();
        REQUIRE (dia[1] == Approx(0.35559999999999997).epsilon (tolerance));
        // check diameter
        auto status = IO->get_pipe_status ();
        REQUIRE (status[0] == true);

        // check length
        auto length = IO->get_pipe_length ();
        REQUIRE (length[1] == Approx(1609.344).epsilon (tolerance));



    }





}