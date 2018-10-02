// Constructor with id and node pointers
pipenetwork::Pipe::Pipe(
    unsigned id, const std::array<std::shared_ptr<pipenetwork::Node>, 2>& nodes)
    : id_{id}, nodes_{nodes} {
  length_ = (nodes_.at(0)->coordinates() - nodes_.at(1)->coordinates()).norm();
}

// Calculate and return discharge using Darcy-Weisbach equation
double pipenetwork::Pipe::discharge_dw() {
  // The Darcy-Weiabach equation is only applicable when heads at both ends of
  // the pipe are known, thus check it
  if (nodes_.at(0)->ishead() && nodes_.at(1)->ishead()) {
    const double dhead = nodes_.at(0)->head() - nodes_.at(1)->head();
    discharge_ = sqrt(std::abs(dhead) * pow(M_PI, 2) * pipenetwork::Gravity(2) *
                      pow(2 * radius_, 5) / (8 * length_ * darcy_friction_));
    // defined flow direction from nodes_.at(0) to nodes.at(1) as positive
    if (dhead < 0) discharge_ *= -1.;
  } else {
    throw std::runtime_error(
        "Unknown head exists, cannot calculate discharge using Darcy Weisbach "
        "equation");
  }
  return discharge_;
}

// Calculate and return discharge using Hazen-Williams equation
double pipenetwork::Pipe::discharge_hw() {
  // To calculate discharge using Hazen-Williams equation, heads at both ends of
  // the pipe should be known, thus check it
  if (nodes_.at(0)->ishead() && nodes_.at(1)->ishead()) {
    const double dhead = nodes_.at(0)->head() - nodes_.at(1)->head();
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
  return discharge_;
}

//! Calculate and return head loss over the pipe using Darcy-Weisbach equation:
double pipenetwork::Pipe::headloss_dw() {
  headloss_ = 8 * length_ * darcy_friction_ * pow(discharge_, 2) /
              (pow(M_PI, 2) * pipenetwork::Gravity(2) * pow(2 * radius_, 5));
  if (discharge_ < 0) headloss_ *= -1.;
  return headloss_;
}

//! Calculate and return headloss over the pipe using Hazen-Williams equation:
double pipenetwork::Pipe::headloss_hw() {
  headloss_ = (10.67 * length_ * pow(discharge_, 1.852)) /
              (pow(pipe_roughness_, 1.852) * pow(2 * radius_, 4.8704));
  if (discharge_ < 0) headloss_ *= -1.;
  return headloss_;
}

//! Return an array of pointers point to the nodes at pipe end
const std::array<std::shared_ptr<const pipenetwork::Node>, 2>
    pipenetwork::Pipe::nodes() {
  std::array<std::shared_ptr<const pipenetwork::Node>, 2> nodes;
  nodes.at(0) = nodes_.at(0);
  nodes.at(1) = nodes_.at(1);
  return nodes;
}
