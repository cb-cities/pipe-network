#include "catch.hpp"
#include <iostream>

#include "node.h"

// Check node class
TEST_CASE("Node is checked", "[Node]") {
  // Tolerance
  const double tolerance = 1.e-6;

  // Index
  unsigned id = 123;
  // Coordinates
  Eigen::Vector3d coords{11.2, 22.3, 33.4};

  // Create a node
  auto node = std::make_unique<pipenetwork::Node>(id, coords);

  // Check node id
  REQUIRE(node->id() == id);

  // Check nodal coordinates
  for (int i = 0; i < coords.size(); ++i)
    REQUIRE(node->coordinates()(i) == Approx(coords(i)).epsilon(tolerance));

  // Check status of discharge before it is initialized
  REQUIRE(node->isdischarge() == false);

  // Check hydraulic head at the node
  SECTION("check hydraulic head at the node") {
    // check PDD polynomial approximation
    double elevation = 110.1;
    node->elevation(elevation);
    auto pdd_coeff1 = node->get_pdd_poly_coef_1();
    auto pdd_coeff2 = node->get_pdd_poly_coef_2();

    REQUIRE(pdd_coeff1[0] == Approx(-18.75).epsilon(tolerance));
    REQUIRE(pdd_coeff1[1] == Approx(6.25).epsilon(tolerance));
    REQUIRE(pdd_coeff1[2] == Approx(1.00009e-12).epsilon(tolerance));
    REQUIRE(pdd_coeff1[3] == Approx(-3.88578e-17).epsilon(tolerance));

    REQUIRE(pdd_coeff2[0] == Approx(-0.624992).epsilon(tolerance));
    REQUIRE(pdd_coeff2[1] == Approx(37.2492).epsilon(tolerance));
    REQUIRE(pdd_coeff2[2] == Approx(-739.978).epsilon(tolerance));
    REQUIRE(pdd_coeff2[3] == Approx(4900.81).epsilon(tolerance));

    // Check nodal elevation
    REQUIRE(node->elevation() == Approx(elevation).epsilon(tolerance));
    // Check if hydraulic head is assigned at the node
    //    REQUIRE(node->ishead() == true);
  }

  // Check discharge at the node
  SECTION("Check discharge at the node") {
    double discharge = -23.4;
    node->demand(discharge);

    // Check discharge at the node
    REQUIRE(node->demand() == Approx(discharge).epsilon(tolerance));
    // Check if discharge is assigned at the node
    REQUIRE(node->isdischarge() == true);
  }

    // Check leak at the node
    SECTION("Check leak at the node") {
      REQUIRE (node->is_leak () == false);
        double leak_d = 0.05;
        node->leak_diameter (leak_d);

        REQUIRE(node->is_leak () == true);
        REQUIRE (node->leak_area () == Approx(0.0019634954).epsilon(tolerance));

        auto leak_coeff = node->get_leak_poly_coef ();
        REQUIRE(leak_coeff[0] == Approx(-97843485.21421291).epsilon(tolerance));
        REQUIRE(leak_coeff[1] == Approx(16307.24753566882).epsilon(tolerance));
        REQUIRE(leak_coeff[2] == Approx(1.000000082740371e-11).epsilon(tolerance));
        REQUIRE(leak_coeff[3] == Approx(1.3791051633906438e-20).epsilon(tolerance));

        node->leak_diameter (0);

        REQUIRE(node->is_leak () == false);


    }

}
