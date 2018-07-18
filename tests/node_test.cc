#include "catch.hpp"

#include "node.h"

//! \brief Check node class
TEST_CASE("Node is checked", "[Node]") {
  // Tolerance
  const double tolerance = 1.e-12;

  // Index
  unsigned id = 123;
  // Coordinates
  Eigen::Vector3d coords{11.2, 22.3, 33.4};

  // Create a node
  auto node = std::make_unique<Node>(id, coords);

  //! Check node id
  REQUIRE(node->id() == id);
  
  //! Check nodal coordinates
  for (int i = 0; i < coords.size(); ++i)
    REQUIRE(node->coordinates()(i) == Approx(coords(i)).epsilon(tolerance));

  //! Check status of head before it is initialized
  REQUIRE(node->ishead() == false);

  //! Check status of discharge before it is initialized
  REQUIRE(node->isdischarge() == false);

  //! Check hydraulic head at the node
  SECTION("check hydraulic head at the node") {
    double head = 110.1;
    node->head(head);

    //! Check hydraulic head at the node
    REQUIRE(node->head() == Approx(head).epsilon(tolerance));
    //! Check if hydraulic head is assigned at the node
    REQUIRE(node->ishead() == true);
  }

  //! Check discharge at the node
  SECTION("Check discharge at the node") {
    double discharge = -23.4;
    node->discharge(discharge);

    //! Check discharge at the node
    REQUIRE(node->discharge() == Approx(discharge).epsilon(tolerance));
    //! Check if discharge is assigned at the node
    REQUIRE(node->isdischarge() == true);
  }

  // TODO: test return nconnections
}

