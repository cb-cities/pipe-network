// Constructor with id, node pointers, diameter, status and max allowable
// velocity
pipenetwork::Pipe::Pipe(
    unsigned id, const std::array<std::shared_ptr<pipenetwork::Node>, 2>& nodes,
    double diameter, bool status, double max_velocity)
    : id_{id},
      nodes_{nodes},
      radius_{diameter / 2.},
      isopen_{status},
      max_velocity_{max_velocity} {
  length_ = (nodes_.at(0)->coordinates() - nodes_.at(1)->coordinates()).norm();
}

// Calculate discharge using Darcy-Weisbach equation
// dhead/length =
// (8*darcy_factor*pow(discharge,2)/(pow(M_PI,2)*g*pow(2*radius,5)) That is,
// discharge =
// sqrt(dhead/length*pow(M_PI,2)*g*pow(2*radius,5)/(8*darcy_friction));
void pipenetwork::Pipe::compute_discharge_darcy_weisbach() {
  // The Darcy-Weiabach equation is only applicable when heads at both ends of
  // the pipe are known, thus check it
  if (nodes_.at(0)->ishead() && nodes_.at(1)->ishead()) {
    double dhead = nodes_.at(0)->head() - nodes_.at(1)->head();
    discharge_ = sqrt(std::abs(dhead) * pow(M_PI, 2) * pipenetwork::Gravity(2) *
                      pow(2 * radius_, 5) / (8. * length_ * darcy_friction_));
    // defined flow direction from nodes_.at(0) to nodes.at(1) as positive
    if (dhead < 0) discharge_ *= -1.;
  } else {
    throw std::runtime_error(
        "Unknown head exists, cannot calculate discharge using Darcy Weisbach "
        "equation");
  }
}

// Calculate discharge using Hazen-Williams equation
// dhead/length =
// (10.67*pow(discharge,1.852)/(pow(pipe_roughness,1.852)*pow(2*radius,4.8704))
// That is, discharge =
// pow((dhead*pow(pipe_roughness,1.852)*pow(2*radius,4.8704)/(10.67*length)),1/1.852);
void pipenetwork::Pipe::compute_discharge_hazen_williams() {
  // To calculate discharge using Hazen-Williams equation, heads at both ends of
  // the pipe should be known, thus check it
  if (nodes_.at(0)->ishead() && nodes_.at(1)->ishead()) {
    double dhead = nodes_.at(0)->head() - nodes_.at(1)->head();
    discharge_ = pow((std::abs(dhead) * pow(pipe_roughness_, 1.852) *
                      pow(2 * radius_, 4.8704) / (10.67 * length_)),
                     1 / 1.852);
    // defined flow direction from nodes_.at(0) to nodes.at(1) as positive
    if (dhead < 0) discharge_ *= -1.;
  } else {
    throw std::runtime_error(
        "Unknown head exists, cannot calculate discharge using Darcy Weisbach "
        "equation");
  }
}

// Calculate head loss over the pipe using Darcy-Weisbach equation:
// dhead =
// length*(8*darcy_factor*pow(discharge,2)/(pow(M_PI,2)*g*pow(2*radius,5))
void pipenetwork::Pipe::compute_headloss_darcy_weisbach() {
  headloss_ = 8. * length_ * darcy_friction_ * pow(discharge_, 2) /
              (pow(M_PI, 2) * pipenetwork::Gravity(2) * pow(2 * radius_, 5));
  if (discharge_ < 0) headloss_ *= -1.;
}

// Calculate headloss over the pipe using Hazen-Williams equation:
// dhead
// = 10.67*length*pow(discharge,1.852)/(pow(pipe_roughness,1.852)*pow(2*radius,4.8704))
void pipenetwork::Pipe::compute_headloss_hazen_williams() {
  headloss_ = (10.67 * length_ * pow(discharge_, 1.852)) /
              (pow(pipe_roughness_, 1.852) * pow(2 * radius_, 4.8704));
  if (discharge_ < 0) headloss_ *= -1.;
}

// Return an array of pointers point to the nodes at pipe end
const std::array<std::shared_ptr<const pipenetwork::Node>, 2>
    pipenetwork::Pipe::nodes() {
  std::array<std::shared_ptr<const pipenetwork::Node>, 2> nodes;
  nodes.at(0) = nodes_.at(0);
  nodes.at(1) = nodes_.at(1);
  return nodes;
}
