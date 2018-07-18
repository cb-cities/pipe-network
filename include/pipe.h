#ifndef PIPE_NETWORK_PIPE_H_
#define PIPE_NETWORK_PIPE_H_

#include "node.h"
#include <array>
#include <math.h>
#include <memory>

//! Pipe class
//! \brief Class that stores the information about pipes
class Pipe {

 public:
  // Constructor with id and node pointers
  //! \param[in] id node id
  //! \param[in] nodes array of node pointers
  Pipe(unsigned id, const std::array<std::shared_ptr<Node>, 2>& nodes)
      : id_{id}, array_node_ptr_{nodes} {}

  //! Destructor
  ~Pipe() { array_node_ptr_.fill(nullptr); }

  //! Copy constructor
  Pipe(const Pipe&) = delete;

  //! Assignment operator
  Pipe& operator=(const Pipe&) = delete;

  //! Move constructor
  Pipe(Pipe&&) = delete;

  //! Return id
  //! \retval id_ id of the pipe
  unsigned id() const { return id_; }

  //! Return array of node pointers
  //! \retval array_node_ptr_ array of node pointers
  std::array<std::shared_ptr<Node>, 2> array_node_ptr() const {
    return array_node_ptr_;
  }

  //! Return discharge in the pipe
  //! \retval discharge_ discharge
  unsigned discharge() const { return discharge_; }

  //! Assign radius of the pipe
  //! \param[in] radius radius of the pipe
  void radius(double radius) { radius_ = radius; }

  //! Assign length of the pipe
  //! \param[in] length of the pipe
  // void length(double length) { length_ = length; }

  //! Assign friction coefficient
  //! \param[in] friction coefficient of the pipe
  // void friction_coef(double friction_coef) { friction_coef_ = friction_coef;
  // }

  //! Assign maximum allowable velocity
  //! \param[in] max_flowrate maximum allowable velocity of flow in the pipe
  void max_velocity(double max_velocity) { max_velocity_ = max_velocity; }

  //! calculate maximum allowable discharge based on ridius and maximum
  //! allowable velocity retval max_discharge maximun allowable discharge
  double max_discharge() {
    double max_discharge = max_velocity_ * M_PI * pow(radius_, 2);
    return max_discharge;
  }

  //! Return pipe broken status
  //! \retval isbroken_ pipe broken status
  bool isbroken() const { return isbroken_; }

 private:
  //! pipe id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! array of node pointers which form the pipe
  std::array<std::shared_ptr<Node>, 2> array_node_ptr_;
  //! discharge in the pipe
  double discharge_{std::numeric_limits<double>::max()};
  //! radius of the pipe
  double radius_{std::numeric_limits<double>::max()};
  //! length of the pipe
  // double length_{std::numeric_limits<double>::max()};
  //! friction coefficient of the pipe
  // double friction_coef_{std::numeric_limits<double>::max()};
  //! maximum allowable velocity of flow in the pipe
  double max_velocity_{std::numeric_limits<double>::max()};
  //! whether the pipe is broken
  bool isbroken_{false};
};

#endif  // PIPE_NETWORK_PIPE_H_
