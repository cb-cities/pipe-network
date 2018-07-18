// Pipe test
#include "node.h"
#include "pipe.h"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

// Check pipe class
TEST_CASE("Pipe is checked", "[Pipe]") {
  // Tolerance
  const double tolerance = 1.e-12;

  // Node index
  unsigned nodeid1 = 100;
  unsigned nodeid2 = 101;

  // Node coordinates
  Eigen::Vector3d coords1(10.0, 20.0, 30.0);
  Eigen::Vector3d coords2(11.1, 21.1, 31.1);

  // Creat an array of two node pointers with previous defined index and
  // coordinates of the nodes
  std::array<std::shared_ptr<Node>, 2> nodeptr;
  nodeptr[0] = std::make_shared<Node>(nodeid1, coords1);
  nodeptr[1] = std::make_shared<Node>(nodeid2, coords2);

  // Pipe index
  unsigned pipeid = 200;

  // Creat a pipe based on previous created node pointers
  auto pipe(std::make_unique<Pipe>(pipeid, nodeptr));

  // Pipe length
  double length = sqrt(3 * pow(1.1, 2));

  // Check pipe id
  REQUIRE(pipe->id() == pipeid);

  // Check pipe length
  REQUIRE(pipe->length() == Approx(length).epsilon(tolerance));

  // Check pipe broken status and initialized status
  REQUIRE(pipe->isbroken() == false);

  // Check radius and max flow velocity of the pipe
  SECTION("check radius and max flow velocity of the pipe") {
    double radius = 10.0;
    double max_velocity = 100.0;
    double discharge = M_PI * 1.e4;
    pipe->radius(radius);
    pipe->max_velocity(max_velocity);

    REQUIRE(pipe->max_discharge() == Approx(discharge).epsilon(tolerance));
  }
}
