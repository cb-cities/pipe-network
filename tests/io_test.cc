#include "catch.hpp"

#include "io.h"

// Check IO class
TEST_CASE("IO is checked", "[IO]") {
  // Tolerance
  const double tolerance = 1.e-12;

  // Create a IO class object
  auto IO = std::make_unique<pipenetwork::IO>();

  // Read node and pipe information from CSV files
  bool read_network =
      IO->read_network("../benchmarks/todini_network_node_iotest.csv",
                       "../benchmarks/todini_network_pipe_iotest.csv");

  // Check read nodal information
  SECTION("Check read nodal information") {
    std::vector<Eigen::Vector3d> nodal_coords = IO->nodal_coordinates();
    std::vector<std::pair<Index, double>> init_nodal_head =
        IO->initial_nodal_head();
    std::vector<std::pair<Index, double>> init_nodal_discharge =
        IO->initial_nodal_discharge();

    // Check coordinates
    REQUIRE(nodal_coords.at(0)(0) == Approx(1.0).epsilon(tolerance));
    REQUIRE(nodal_coords.at(0)(1) == Approx(3.0).epsilon(tolerance));
    REQUIRE(nodal_coords.at(0)(2) == Approx(0.0).epsilon(tolerance));
    REQUIRE(nodal_coords.at(1)(0) == Approx(0.0).epsilon(tolerance));
    REQUIRE(nodal_coords.at(1)(1) == Approx(2.0).epsilon(tolerance));
    REQUIRE(nodal_coords.at(1)(2) == Approx(0.0).epsilon(tolerance));
    REQUIRE(nodal_coords.at(2)(0) == Approx(2.0).epsilon(tolerance));
    REQUIRE(nodal_coords.at(2)(1) == Approx(2.0).epsilon(tolerance));
    REQUIRE(nodal_coords.at(2)(2) == Approx(0.0).epsilon(tolerance));
    REQUIRE(nodal_coords.at(3)(0) == Approx(0.0).epsilon(tolerance));
    REQUIRE(nodal_coords.at(3)(1) == Approx(0.0).epsilon(tolerance));
    REQUIRE(nodal_coords.at(3)(2) == Approx(0.0).epsilon(tolerance));
    REQUIRE(nodal_coords.at(4)(0) == Approx(2.0).epsilon(tolerance));
    REQUIRE(nodal_coords.at(4)(1) == Approx(0.0).epsilon(tolerance));
    REQUIRE(nodal_coords.at(4)(2) == Approx(0.0).epsilon(tolerance));

    // Check initial nodal head
    REQUIRE(init_nodal_head.at(0).first == 0);
    REQUIRE(init_nodal_head.at(0).second == Approx(100.0).epsilon(tolerance));
    REQUIRE(init_nodal_head.at(1).first == 1);
    REQUIRE(init_nodal_head.at(1).second == Approx(99.0).epsilon(tolerance));
    REQUIRE(init_nodal_head.at(2).first == 2);
    REQUIRE(init_nodal_head.at(2).second == Approx(98.0).epsilon(tolerance));
    REQUIRE(init_nodal_head.at(3).first == 3);
    REQUIRE(init_nodal_head.at(3).second == Approx(97.0).epsilon(tolerance));

    // Check initial nodal discharge
    REQUIRE(init_nodal_discharge.at(0).first == 0);
    REQUIRE(init_nodal_discharge.at(0).second ==
            Approx(-100.0).epsilon(tolerance));
    REQUIRE(init_nodal_discharge.at(1).first == 1);
    REQUIRE(init_nodal_discharge.at(1).second ==
            Approx(10.0).epsilon(tolerance));
    REQUIRE(init_nodal_discharge.at(2).first == 2);
    REQUIRE(init_nodal_discharge.at(2).second ==
            Approx(20.0).epsilon(tolerance));
    REQUIRE(init_nodal_discharge.at(3).first == 3);
    REQUIRE(init_nodal_discharge.at(3).second ==
            Approx(30.0).epsilon(tolerance));
    REQUIRE(init_nodal_discharge.at(4).first == 4);
    REQUIRE(init_nodal_discharge.at(4).second ==
            Approx(40.0).epsilon(tolerance));
  }

  // Check read pipe information
  SECTION("Check read pipe information") {
    std::vector<std::pair<Index, Index>> node_pairs = IO->node_pairs();
    std::vector<double> diameters = IO->diameters();
    std::vector<double> roughness = IO->roughness();
    std::vector<bool> pipe_status = IO->pipe_status();
    std::vector<std::pair<Index, double>> init_pipe_discharge =
        IO->initial_pipe_discharge();

    // Check start and end nodes
    REQUIRE(node_pairs.at(0).first == 0);
    REQUIRE(node_pairs.at(0).second == 1);
    REQUIRE(node_pairs.at(1).first == 0);
    REQUIRE(node_pairs.at(1).second == 2);
    REQUIRE(node_pairs.at(2).first == 1);
    REQUIRE(node_pairs.at(2).second == 2);
    REQUIRE(node_pairs.at(3).first == 1);
    REQUIRE(node_pairs.at(3).second == 3);
    REQUIRE(node_pairs.at(4).first == 1);
    REQUIRE(node_pairs.at(4).second == 4);
    REQUIRE(node_pairs.at(5).first == 2);
    REQUIRE(node_pairs.at(5).second == 4);
    REQUIRE(node_pairs.at(6).first == 3);
    REQUIRE(node_pairs.at(6).second == 4);

    // Check pipe diameters
    REQUIRE(diameters.at(0) == Approx(0.1).epsilon(tolerance));
    REQUIRE(diameters.at(1) == Approx(0.2).epsilon(tolerance));
    REQUIRE(diameters.at(2) == Approx(0.3).epsilon(tolerance));
    REQUIRE(diameters.at(3) == Approx(0.4).epsilon(tolerance));
    REQUIRE(diameters.at(4) == Approx(0.5).epsilon(tolerance));
    REQUIRE(diameters.at(5) == Approx(0.6).epsilon(tolerance));
    REQUIRE(diameters.at(6) == Approx(0.7).epsilon(tolerance));

    // Check pipe roughness
    REQUIRE(roughness.at(0) == Approx(147683.3637).epsilon(tolerance));
    REQUIRE(roughness.at(1) == Approx(4102.777965).epsilon(tolerance));
    REQUIRE(roughness.at(2) == Approx(1238.160737).epsilon(tolerance));
    REQUIRE(roughness.at(3) == Approx(1598.561076).epsilon(tolerance));
    REQUIRE(roughness.at(4) == Approx(430.5609768).epsilon(tolerance));
    REQUIRE(roughness.at(5) == Approx(137.5889198).epsilon(tolerance));
    REQUIRE(roughness.at(6) == Approx(133.3735196).epsilon(tolerance));

    // Check pipe status
    for (int i = 0; i < pipe_status.size(); i++)
      REQUIRE(pipe_status.at(i) == true);

    // Check initial pipe discharge
    REQUIRE(init_pipe_discharge.at(0).first == 0);
    REQUIRE(init_pipe_discharge.at(0).second ==
            Approx(80.0).epsilon(tolerance));
    REQUIRE(init_pipe_discharge.at(1).first == 1);
    REQUIRE(init_pipe_discharge.at(1).second ==
            Approx(20.0).epsilon(tolerance));
    REQUIRE(init_pipe_discharge.at(2).first == 2);
    REQUIRE(init_pipe_discharge.at(2).second ==
            Approx(10.0).epsilon(tolerance));
    REQUIRE(init_pipe_discharge.at(3).first == 3);
    REQUIRE(init_pipe_discharge.at(3).second ==
            Approx(40.0).epsilon(tolerance));
    REQUIRE(init_pipe_discharge.at(4).first == 4);
    REQUIRE(init_pipe_discharge.at(4).second ==
            Approx(20.0).epsilon(tolerance));
    REQUIRE(init_pipe_discharge.at(5).first == 6);
    REQUIRE(init_pipe_discharge.at(5).second ==
            Approx(10.0).epsilon(tolerance));
  }
}
