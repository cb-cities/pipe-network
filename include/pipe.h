#ifndef PIPE_NETWORK_PIPE_H_
#define PIPE_NETWORK_PIPE_H_

#include <cmath>

#include <array>
#include <exception>
#include <memory>

#include "node.h"

//! Pipe class
//! \brief Class that stores the information about pipes
class Pipe {

 public:
  // Constructor with id and node pointers
  //! \param[in] id node id
  //! \param[in] nodes array of node pointers
  Pipe(unsigned id, const std::array<std::shared_ptr<Node>, 2>& nodes);

  //! Destructor
  ~Pipe() { nodes_.fill(nullptr); }

  //! Copy constructor
  Pipe(const Pipe&) = delete;

  //! Assignment operator
  Pipe& operator=(const Pipe&) = delete;

  //! Move constructor
  Pipe(Pipe&&) = delete;

  //! Return id
  //! \retval id_ id of the pipe
  unsigned id() const { return id_; }

  //! Assign radius of the pipe
  //! \param[in] radius radius of the pipe
  void radius(double radius) { radius_ = radius; }

  //! Return length of the pipe
  //! retval length_ pipe length
  double length() { return length_; }

  //! Assign Dracy friction factor
  //! \param[in] darcy_friction Darcy friction factor of the pipe
  void darcy_friction(double darcy_friction) {
    darcy_friction_ = darcy_friction;
  }

  //! Calculate and return discharge using Darcy-Weisbach equation:
  //! dhead = (8*darcy_factor*pow(discharge,2)/(pow(M_PI,2)*g*pow(2*radius,5))
  //! That is, discharge =
  //! sqrt(dhead*pow(M_PI,2)*g*pow(2*radius,5)/(8*darcy_friction)); \retval
  //! discharge_ discharge in the pipe
  double discharge();

  //! Assign maximum allowable velocity
  //! \param[in] max_flowrate maximum allowable velocity of flow in the pipe
  void max_velocity(double max_velocity) { max_velocity_ = max_velocity; }

  //! calculate maximum allowable discharge based on ridius and maximum
  //! allowable velocity \retval max_discharge maximun allowable discharge
  double max_discharge() { return (max_velocity_ * M_PI * pow(radius_, 2)); }

  //! Return pipe broken status
  //! \retval isbroken_ pipe broken status
  bool isbroken() const { return isbroken_; }

 private:
  //! pipe id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! array of node pointers which form the pipe
  std::array<std::shared_ptr<Node>, 2> nodes_;
  //! discharge in the pipe
  double discharge_{std::numeric_limits<double>::max()};
  //! radius of the pipe
  double radius_{std::numeric_limits<double>::max()};
  //! length of the pipe
  double length_;
  //! Darcy friction factor of the pipe
  double darcy_friction_{std::numeric_limits<double>::max()};
  //! maximum allowable velocity of flow in the pipe
  double max_velocity_{std::numeric_limits<double>::max()};
  //! whether the pipe is broken
  bool isbroken_{false};
  //! gravitational acceleration
  const double g_{9.81};
};

#endif  // PIPE_NETWORK_PIPE_H_
