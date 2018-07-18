#ifndef PIPE_NETWORK_PIPE_H_
#define PIPE_NETWORK_PIPE_H_

#include <cmath>

#include <array>
#include <memory>

#include "node.h"

//! Pipe class
//! \brief Class that stores the information about pipes
class Pipe {

 public:
  // Constructor with id and node pointers
  //! \param[in] id node id
  //! \param[in] nodes array of node pointers
  Pipe(unsigned id, const std::array<std::shared_ptr<Node>, 2>& nodes)
      : id_{id}, nodes_{nodes} {
    Eigen::Vector3d distance =
        nodes_.at(0)->coordinates() - nodes_.at(1)->coordinates();
    length_ = distance.norm();
  }

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

  //! Return discharge in the pipe
  //! \retval discharge_ discharge
  unsigned discharge() const { return discharge_; }

  //! Assign radius of the pipe
  //! \param[in] radius radius of the pipe
  void radius(double radius) { radius_ = radius; }

  //! Return length of the pipe
  //! retval length_ pipe length
  double length() { return length_; }

  //! Assign friction coefficient
  //! \param[in] friction coefficient of the pipe
  // void friction_coef(double friction_coef) { friction_coef_ = friction_coef;
  // }

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
  //! friction coefficient of the pipe
  // double friction_coef_{std::numeric_limits<double>::max()};
  //! maximum allowable velocity of flow in the pipe
  double max_velocity_{std::numeric_limits<double>::max()};
  //! whether the pipe is broken
  bool isbroken_{false};
};

#endif  // PIPE_NETWORK_PIPE_H_
