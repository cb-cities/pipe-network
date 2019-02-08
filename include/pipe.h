#ifndef PIPE_NETWORK_PIPE_H_
#define PIPE_NETWORK_PIPE_H_

#include <cmath>

#include <array>
#include <exception>
#include <memory>

#include "node.h"
#include "settings.h"

namespace pipenetwork {

//! Pipe class
//! \brief Class that stores the information about pipes
class Pipe {

 public:
  //! Constructor with pipe id, node pointers, diameter, status and max
  //! allowable velocity
  //! \param[in] id pipe id
  //! \param[in] nodes array of node pointers
  //! \param[in] diameter of the pipe
  //! \param[in] roughness pipe roughness for Hazen-Williams equation
  //! \param[in] status of the pipe, true for open, flase for close
  //! \param[in] maximum allowable velocity
  Pipe(unsigned id,
       const std::array<std::shared_ptr<pipenetwork::Node>, 2>& nodes,
       double diameter, double roughness, bool status,
       double max_velocity = std::numeric_limits<double>::max());

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

  //! Return length of the pipe
  //! \retval length_ pipe length
  double length() { return length_; }

  //! Assign Dracy friction factor
  //! \param[in] darcy_friction Darcy friction factor of the pipe
  void darcy_friction(double darcy_friction) {
    darcy_friction_ = darcy_friction;
  }

  //! Assign discharge to the pipe
  //! \param[in] iter_discharge dischargeof the pipe during iteration
  void iter_discharge(double iter_discharge) {
    iter_discharge_ = iter_discharge;
  }

  //! Initialize discharge with input value
  //! \param[in] discharge input discharge value of the pipe
  void initialize_discharge(double discharge = 0.001) {
    discharge_ = discharge;
    iter_discharge_ = discharge_;
  }

  //! Calculate discharge using Darcy-Weisbach equation
  //! Calculate discharge in m^3/s from head difference in m
  //! SI unit meter and second are used in the whole equation
  void compute_discharge_darcy_weisbach();

  //! Calculate discharge using Hazen-Williams equation
  //! Calculate discharge in m^3/s from head difference in m
  //! SI unit meter and second are used in the whole equation
  void compute_discharge_hazen_williams();

  //! Calculate and return derivative of Hazen-Williams equation with respect to
  //! pipe discharge SI unit meter and second are used in the whole equation
  //! \retval derivative of Hazen-Williams equation with respect to pipe
  //! discharge
  double deriv_hazen_williams_discharge();

  //! Return calculated discharge in the pipe
  //! \retval discharge_ discharge in the pipe
  double discharge() const { return discharge_; }

  //! Calculate head loss over the pipe using Darcy-Weisbach equation:
  //! Calculate head loss in m from discharge in m^3/s
  //! SI unit meter and second are used in the whole equation
  void compute_headloss_darcy_weisbach();

  //! Calculate headloss over the pipe using Hazen-Williams equation:
  //! Calculate head loss in m from discharge in m^3/s
  //! SI unit meter and second are used in the whole equation
  void compute_headloss_hazen_williams();

  //! Return calculated headloss over the pipe
  //! \retval headloss_ headloss over the pipe
  double headloss() const { return headloss_; }

  //! Return assigned discharge during iteration of the pipe
  //! \retval iter_discharge_ discharge of the pipe during iteration
  double iter_discharge() const { return iter_discharge_; }

  //! calculate maximum allowable discharge based on ridius and maximum
  //! allowable velocity \retval max_discharge maximun allowable discharge
  double max_discharge() { return (max_velocity_ * M_PI * pow(radius_, 2)); }

  //! Return pipe open status
  //! \retval isopen_ whether pipe open or close
  bool isopen() const { return isopen_; }

  //! Return pipe broken status
  //! \retval isbroken_ pipe broken status
  bool isbroken() const { return isbroken_; }

  //! Return an array of pointers point to the nodes at pipe end
  const std::array<std::shared_ptr<const pipenetwork::Node>, 2> nodes();

 private:
  //! pipe id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! array of node pointers which form the pipe
  std::array<std::shared_ptr<pipenetwork::Node>, 2> nodes_;
  //! discharge in the pipe
  double discharge_{std::numeric_limits<double>::max()};
  //! discharge assigned to the pipe during iteration
  double iter_discharge_{std::numeric_limits<double>::max()};
  //! headloss over the pipe
  double headloss_{std::numeric_limits<double>::max()};
  //! radius of the pipe
  double radius_{std::numeric_limits<double>::max()};
  //! length of the pipe
  double length_;
  //! Darcy friction factor of the pipe used in Darcy-Weisbach equation
  double darcy_friction_{std::numeric_limits<double>::max()};
  //! pipe roughness coefficient used in Hazen-Williams equation
  double pipe_roughness_{std::numeric_limits<double>::max()};
  //! maximum allowable velocity of flow in the pipe
  double max_velocity_{std::numeric_limits<double>::max()};
  //! whether the pipe is broken
  bool isbroken_{false};
  //! whether the pipe is open
  bool isopen_{true};
};
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_PIPE_H_
