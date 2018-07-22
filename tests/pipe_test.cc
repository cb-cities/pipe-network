#include "catch.hpp"

#include "node.h"
#include "pipe.h"

// Check pipe class
TEST_CASE("Pipe is checked", "[Pipe]") {
  // Tolerance
  const double tolerance = 1.e-12;

  // Node indices
  const unsigned nodeid1 = 100;
  const unsigned nodeid2 = 101;

  // Nodal coordinates
  const Eigen::Vector3d coords1(10.0, 20.0, 30.0);
  const Eigen::Vector3d coords2(11.1, 21.1, 31.1);

  // Creat an array of two node pointers with previous defined index and
  // coordinates of the nodes
  std::array<std::shared_ptr<Node>, 2> nodes;
  nodes[0] = std::make_shared<Node>(nodeid1, coords1);
  nodes[1] = std::make_shared<Node>(nodeid2, coords2);

  // Pipe index
  const unsigned pipeid = 200;

  // Creat a pipe based on previous created node pointers
  auto pipe = std::make_unique<Pipe>(pipeid, nodes);

  // Pipe length
  const double length = sqrt(3 * pow(1.1, 2));

  // Check pipe id
  REQUIRE(pipe->id() == pipeid);

  // Check pipe length
  REQUIRE(pipe->length() == Approx(length).epsilon(tolerance));

  // Check pipe broken status and initialized status
  REQUIRE(pipe->isbroken() == false);

  // Ckeck node in pipe
  REQUIRE(pipe->isnode(nodes[0]) == true);
  REQUIRE(pipe->isnode(nodes[1]) == true);

  // Check coordinates at pipe ends
  REQUIRE(pipe->end_coordinates().at(0) == nodes.at(0)->coordinates());
  REQUIRE(pipe->end_coordinates().at(1) == nodes.at(1)->coordinates());

  // Check radius, discharge, max flow velocity and Darcy friction factor of the
  // pipe
  SECTION(
      "Check radius, discharge, max flow velocity and Darcy friction factor of "
      "the pipe") {
    const double radius = 10.0;
    const double max_velocity = 100.0;
    const double max_discharge = M_PI * 1.e4;
    const double head1 = 110.0;
    const double head2 = 100.0;
    const double darcy_friction = 0.1;
    const double discharge =
        sqrt(10 * pow(M_PI, 2) * 9.81 * pow(2 * 10, 5) / (8 * 0.1));
    nodes[0]->head(head1);
    nodes[1]->head(head2);
    pipe->radius(radius);
    pipe->max_velocity(max_velocity);
    pipe->darcy_friction(darcy_friction);

    // Check radius and max flow velocity of the pipe
    REQUIRE(pipe->max_discharge() == Approx(max_discharge).epsilon(tolerance));
    // Check discharge and Darcy friction factor of the pipe
    REQUIRE(pipe->discharge() == Approx(discharge).epsilon(tolerance));
  }
}
