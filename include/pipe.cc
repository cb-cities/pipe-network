// Constructor with id and node pointers
//! \param[in] id node id
//! \param[in] nodes array of node pointers
Pipe::Pipe(unsigned id, const std::array<std::shared_ptr<Node>, 2>& nodes)
    : id_{id}, nodes_{nodes} {
  Eigen::Vector3d distance =
      nodes_.at(0)->coordinates() - nodes_.at(1)->coordinates();
  length_ = distance.norm();
}

//! Calculate and return discharge using Darcy-Weisbach equation
//! retval discharge_ discharge in the pipe
double Pipe::discharge() {
  if (nodes_[0]->ishead() && nodes_[1]->ishead()) {
    double dhead = abs(nodes_[0]->head() - nodes_[1]->head());
    discharge_ = sqrt(dhead * pow(M_PI, 2) * g_ * pow(2 * radius_, 5) /
                      (8 * darcy_friction_));
    return discharge_;
  } else {
    throw std::runtime_error(
        "Unknown head exists, cannot calculate discharge using Darcy Weisbach "
        "equation");
  }
}
