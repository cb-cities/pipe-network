#include "catch.hpp"

#include "node.h"
#include "pipe.h"
#include "settings.h"

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
  std::array<std::shared_ptr<pipenetwork::Node>, 2> nodes;
  nodes.at(0) = std::make_shared<pipenetwork::Node>(nodeid1, coords1);
  nodes.at(1) = std::make_shared<pipenetwork::Node>(nodeid2, coords2);

  // Pipe index
  const unsigned pipeid = 201;

  // Dismeter of the pipe in m
  const double diameter = 20.0;

  // Pipe open status
  bool status = true;

  // Creat pipes based on previous created node pointers
  auto pipe1 =
      std::make_unique<pipenetwork::Pipe>(pipeid, nodes, diameter, status);

  // Pipe length
  // To calculate length between two points in 3d space
  // length = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2))
  const double length = sqrt(3 * pow(1.1, 2));

  // Check pipe id
  REQUIRE(pipe1->id() == pipeid);

  // Check pipe length
  REQUIRE(pipe1->length() == Approx(length).epsilon(tolerance));

  // Check initialized discharge and return discharge
  pipe1->initialize_discharge();
  REQUIRE(pipe1->discharge() == Approx(0.001).epsilon(tolerance));

  // Check pipe open status
  REQUIRE(pipe1->isopen() == true);

  // Check pipe broken status and initialized status
  REQUIRE(pipe1->isbroken() == false);

  // Check radius, discharge, headloss, max flow velocity, pipe roughness and
  // Darcy friction factor of the pipe
  SECTION(
      "Check radius, discharge, max flow velocity and Darcy friction factor of "
      "the pipe") {

    // Maximum allowable flow velocity of the pipe in m/s
    const double max_velocity = 100.0;

    // Creat pipes based on previous created node pointers
    auto pipe2 = std::make_unique<pipenetwork::Pipe>(pipeid, nodes, diameter,
                                                     status, max_velocity);

    // Initialize discharge in the pipe and check
    const double init_discharge = 10.0;
    pipe2->initialize_discharge(init_discharge);
    REQUIRE(pipe2->discharge() == Approx(init_discharge).epsilon(tolerance));

    // Maximum allowable discharge of the pipe in m3/s
    // Calculated by max_discharge=M_PI*pow(radius,2)*max_velocity
    const double max_discharge = M_PI * 1.e4;

    // Water head at two nodes in m
    const double head1 = 110.0;
    const double head2 = 100.0;

    // Pipe length
    const double length = sqrt(3 * pow(1.1, 2));

    // Gravity
    const double gravity = 9.81;

    // Dimensionless Darcy friction factor (for Darcy-Weisbach equation)
    const double darcy_friction = 0.1;

    // Pipe roughness coefficient (for Hazen-Williams equation)
    const double pipe_roughness = 100;

    // Assign defined variables to nodes and pipe
    nodes.at(0)->head(head1);
    nodes.at(1)->head(head2);
    pipe2->darcy_friction(darcy_friction);
    pipe2->pipe_roughness(pipe_roughness);

    // Check radius and max flow velocity of the pipe
    REQUIRE(pipe2->max_discharge() == Approx(max_discharge).epsilon(tolerance));

    // Check discharge and Darcy friction factor of the pipe using
    // Darcy-Weisbach equation. Compared with hand calculation:
    // discharge = sqrt{\frac{(110.0-100.0) \times \pi^2 \times 9.81
    // \times 20.0^5}{8.0 \times 0.1 \times squt{3 \times 1.1^2}}}
    // = 45085.585688339
    pipe2->compute_discharge_darcy_weisbach();
    REQUIRE(pipe2->discharge() == Approx(45085.585688339).epsilon(tolerance));
    // Check headloss of the pipe using Darcy-Weisbach equation.
    // Compared with hand calculation:
    // headloss = \frac{8 \times sqrt{3 \times 1.1^2} \times 0.1 \times
    // 45085.585688339^2}{pi^2 \times 9.81 \times 20.0^5} = 10.0
    pipe2->compute_headloss_darcy_weisbach();
    REQUIRE(pipe2->headloss() == Approx(10.0).epsilon(tolerance));

    // Check discharge and pipe roughness coefficient of the pipe using
    // Hazen-Williams equation. Compared with hand calculation:
    // discharge = (\frac{(110.0-100) \times 100^1.852 \times  20.0^4.8704}
    // {10.67 \times sqrt{3 \times 1.1^2}})^{\frac{1}{1.852}}
    // = 179922.60192865
    pipe2->compute_discharge_hazen_williams();
    REQUIRE(pipe2->discharge() == Approx(179922.60192865).epsilon(tolerance));
    // Check headloss of the pipe using Hazen-Williams equation.
    // Compared with hand calculation:
    // Hand calculation: headloss = \frac{10.67 \times sqrt{3 \times 1.1^2}
    // \times 179922.60192865^1.852}{100^1.852 \times 20.0^4.8704} = 10.0
    pipe2->compute_headloss_hazen_williams();
    REQUIRE(pipe2->headloss() == Approx(10.0).epsilon(tolerance));
  }

  // Check return pointer to const Node
  SECTION("Check return pointer to const Node") {

    // Check return ids of the nodes
    REQUIRE(pipe1->nodes().at(0)->id() == 100);
    REQUIRE(pipe1->nodes().at(1)->id() == 101);
    // Check return head assignment status
    REQUIRE(pipe1->nodes().at(0)->ishead() == false);
    REQUIRE(pipe1->nodes().at(1)->ishead() == false);
  }
}
