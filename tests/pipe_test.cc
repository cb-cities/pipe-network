// Pipe test
#include "node.h"
#include "pipe.h"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

//! \brief Ckeck pipe class
TEST_CASE("Pipe is checked", "[Pipe]") {
  unsigned nodeid1 = 100;
  unsigned nodeid2 = 101;
  Eigen::Vector3d coords1(10.0, 20.0, 30.0);
  Eigen::Vector3d coords2(11.1, 21.1, 31.1);
  std::array<std::shared_ptr<Node>, 2> nodeptr;
  nodeptr[0] = std::make_shared<Node>(nodeid1, coords1);
  nodeptr[1] = std::make_shared<Node>(nodeid2, coords2);
  unsigned pipeid = 200;
  auto pipe(std::make_unique<Pipe>(pipeid, nodeptr));
  const double tolerance = 1.e-12;

  //! Check pipe id
  REQUIRE(pipe->id() == pipeid);
  //! Check array of node pointers
  REQUIRE(pipe->array_node_ptr()[0]->id() == nodeid1);
  REQUIRE(pipe->array_node_ptr()[1]->id() == nodeid2);
  //! Check pipe broken status and initialized status
  REQUIRE(pipe->isbroken() == false);

  //! Check radius, discharge and max flow velocity of the pipe
  SECTION("check radius of the pipe") {
    double radius = 10.0;
    double max_velocity = 100.0;
    double discharge = M_PI * 1.e4;
    pipe->radius(radius);
    pipe->max_velocity(max_velocity);

    REQUIRE(pipe->max_discharge() == Approx(discharge).epsilon(tolerance));
  }
}
