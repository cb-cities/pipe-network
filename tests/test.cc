//Node test
#include "node.h"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

//! \brief Ckeck node class
TEST_CASE("Node is checked", "[Node]") {
  const double tolerance = 1.e-12;
  unsigned id = 123;
  Eigen::Vector3d coords(11.2, 22.3, 33.4);
  auto node(std::make_unique<Node>(id, coords));
  
  //! Check function to return node id
  REQUIRE(node->id() == id);
  //! Check function to return node coordinates
  REQUIRE(node->coordinates()(0) == Approx(coords(0)).epsilon(tolerance));
  REQUIRE(node->coordinates()(1) == Approx(coords(1)).epsilon(tolerance));
  REQUIRE(node->coordinates()(2) == Approx(coords(2)).epsilon(tolerance));
  //! Check function to return head assignment status and initialized status
  REQUIRE(node->ishead() == false);
  //! Check function to return discharge assignment status and initialized status
  REQUIRE(node->isdischarge() == false);
  
  //! Check functions to assign and return hydraulic head at the node
  SECTION("check hydraulic head at the node") {
    double head = 110.1;
    node->head(head);
    
    //! Check funtion to return hydraulic head at the node
    REQUIRE(node->head() == Approx(head).epsilon(tolerance));
    //! Check function to assign hydraulic head at the node
    REQUIRE(node->ishead() == true);
  }
  
  //! Check functions to assign and return discharge at the node
  SECTION("Check discharge at the node") {
    double discharge = -23.4;
    node->discharge(discharge);
    
    //! Check function to return discharge at the node
    REQUIRE(node->discharge() == Approx(discharge).epsilon(tolerance));
    //! Check function to assign discharge at the node
    REQUIRE(node->isdischarge() == true);
  }

  // test return nconnections
}
