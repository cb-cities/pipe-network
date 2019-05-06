#include "catch.hpp"

#include "node.h"
#include <iostream>

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

  // TODO: test return nconnections
}
