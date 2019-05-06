#include "pipe.h"

// Constructor with id, node pointers, diameter, status and max allowable
// velocity
pipenetwork::Pipe::Pipe(
    unsigned id, const std::array<std::shared_ptr<pipenetwork::Node>, 2>& nodes,
    double diameter, double roughness, bool status, double max_velocity)
    : id_{id},
      nodes_{nodes},
      radius_{diameter / 2.},
      pipe_roughness_{roughness},
      isopen_{status},
      max_velocity_{max_velocity} {
  length_ = (nodes_.at(0)->coordinates() - nodes_.at(1)->coordinates()).norm();
}

// Calculate discharge using Darcy-Weisbach equation
// $discharge = (\frac{dhead \times \pi^2 \times g \times (2radius)^5}{8 \times
// length \times darcy_friction})^0.5$
// SI unit meter and second are used in the whole equation
void pipenetwork::Pipe::compute_discharge_darcy_weisbach() {
  const double dhead = nodes_.at(0)->head() - nodes_.at(1)->head();
  discharge_ = sqrt(std::abs(dhead) * pow(M_PI, 2) * pipenetwork::Gravity(2) *
                    pow(2 * radius_, 5) / (8. * length_ * darcy_friction_));
  // defined flow direction from nodes_.at(0) to nodes.at(1) as positive
  if (dhead < 0) discharge_ *= -1.;
}

// Calculate discharge using Hazen-Williams equation
// discharge = (\frac{dhead \times pipe_roughness^1.852 \times
// (2radius)^4.8704}{10.67 \times length})^(\frac{1}{1.852})
// SI unit meter and second are used in the whole equation
void pipenetwork::Pipe::compute_discharge_hazen_williams() {
  const double dhead = nodes_.at(0)->head() - nodes_.at(1)->head();
  discharge_ = pow((std::abs(dhead) * pow(pipe_roughness_, 1.852) *
                    pow(2 * radius_, 4.8704) / (10.67 * length_)),
                   1 / 1.852);
  // defined flow direction from nodes_.at(0) to nodes.at(1) as positive
  if (dhead < 0) discharge_ *= -1.;
}

// Calculate and return derivative of Hazen-Williams equation with respect to
// pipe discharge Headloss equation (Hazen-Williams): Residual =
// (Start_node_head - end_node_head) - \frac{10.67 \times length \times
//            pipe_discharge^1.852}{pipe_roughness^1.852 \times
//            (2radius)^4.8704}
// Take derivative with respect to pipe_discharge:
// \frac{-1.852 \times 10.67 \times length}{pow(pipe_roughness,1.852)
// \times pow(2radius,4.8704)} \times pow(pipe_discharge,0.852).
// SI unit meter and second are used in the whole equation
double pipenetwork::Pipe::deriv_hazen_williams_discharge() {
  double coeff = 10.67 * length_ /
                 (pow(pipe_roughness_, 1.852) * pow(2 * radius_, 4.8704));
  double deriv = -1.852 * coeff * pow(std::abs(iter_discharge_), 0.852);
  if (iter_discharge_ < 0) deriv *= -1.;
  return deriv;
}

// Calculate head loss over the pipe using Darcy-Weisbach equation:
// dhead = \frac{8 \times length \times darcy_factor \times discharge^2}{\pi^2
// \times g \times (2radius)^5}
// SI unit meter and second are used in the whole equation
void pipenetwork::Pipe::compute_headloss_darcy_weisbach() {
  headloss_ = 8. * length_ * darcy_friction_ * pow(iter_discharge_, 2) /
              (pow(M_PI, 2) * pipenetwork::Gravity(2) * pow(2 * radius_, 5));
  if (iter_discharge_ < 0) headloss_ *= -1.;
}

// Calculate headloss over the pipe using Hazen-Williams equation:
// dhead = \frac{10.67 \times length \times
// discharge^1.852}{pipe_roughness^1.852 \times (2radius)^4.8704}
// SI unit meter and second are used in the whole equation
void pipenetwork::Pipe::compute_headloss_hazen_williams() {
  headloss_ = (10.67 * length_ * pow(std::abs(iter_discharge_), 1.852)) /
              (pow(pipe_roughness_, 1.852) * pow(2 * radius_, 4.8704));
  if (iter_discharge_ < 0) headloss_ *= -1.;
}

// Return an array of pointers point to the nodes at pipe end
const std::array<std::shared_ptr<const pipenetwork::Node>, 2>
    pipenetwork::Pipe::nodes() {
  std::array<std::shared_ptr<const pipenetwork::Node>, 2> nodes;
  nodes.at(0) = nodes_.at(0);
  nodes.at(1) = nodes_.at(1);
  return nodes;
}
