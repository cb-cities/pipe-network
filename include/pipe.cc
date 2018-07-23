// Constructor with id and node pointers
pipenetwork::Pipe::Pipe(
    unsigned id, const std::array<std::shared_ptr<pipenetwork::Node>, 2>& nodes)
    : id_{id}, nodes_{nodes} {
  length_ = (nodes_.at(0)->coordinates() - nodes_.at(1)->coordinates()).norm();
}

// Calculate and return discharge using Darcy-Weisbach equation
double pipenetwork::Pipe::discharge() {
  // The Darcy-Weiabach equation is only applicable when heads at both ends of
  // the pipe are known, thus check it
  if (nodes_.at(0)->ishead() && nodes_.at(1)->ishead()) {
    const double dhead = nodes_.at(0)->head() - nodes_.at(1)->head();
    discharge_ = sqrt(std::abs(dhead) * pow(M_PI, 2) * pipenetwork::Gravity(2) *
                      pow(2 * radius_, 5) / (8 * darcy_friction_));
    // defined flow direction from nodes_.at(0) to nodes.at(1) as positive
    if (dhead < 0) discharge_ *= -1.;
  } else {
    throw std::runtime_error(
        "Unknown head exists, cannot calculate discharge using Darcy Weisbach "
        "equation");
  }
  return discharge_;
}
