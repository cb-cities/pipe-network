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
  //! Constructor with pipe id, node pointers, diameter and status
  //! \param[in] id pipe id
  //! \param[in] nodes array of node pointers
  //! \param[in] diameter of the pipe
  //! \param[in] status of the pipe, true for open, flase for close
  Pipe(unsigned id,
       const std::array<std::shared_ptr<pipenetwork::Node>, 2>& nodes,
       const double diameter, bool status);

  //! Constructor with pipe id, node pointers, diameter, status and max
  //! allowable velocity \param[in] id pipe id 
  //! \param[in] nodes array of nodepointers 
  //! \param[in] diameter of the pipe 
  //! \param[in] status of the pipe, true for open, flase for close 
  //! \param[in] maximum allowable velocity
  Pipe(unsigned id,
       const std::array<std::shared_ptr<pipenetwork::Node>, 2>& nodes,
       const double diameter, bool status, const double max_velocity);

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

  //! Assign pipe roughness coefficient
  //! \param[in] pipe_roughness roughness coefficient of the pipe
  void pipe_roughness(double pipe_roughness) {
    pipe_roughness_ = pipe_roughness;
  }

  //! Initialize discharge to a small value
  void initialize_discharge() { discharge_ = 0.001; }

  //! Initialize discharge with input value
  //! \param[in] discharge input discharge value of the pipe
  void initialize_discharge(double discharge) { discharge_ = discharge; }

  //! Calculate discharge using Darcy-Weisbach equation
  void compute_discharge_darcy_weisbach();

  //! Calculate discharge using Hazen-Williams equation
  void compute_discharge_hazen_williams();

  //! Return calculated discharge in the pipe
  //! \retval discharge_ discharge in the pipe
  double discharge() const { return discharge_; }

  //! Calculate head loss over the pipe using Darcy-Weisbach equation:
  void compute_headloss_darcy_weisbach();

  //! Calculate headloss over the pipe using Hazen-Williams equation:
  void compute_headloss_hazen_williams();

  //! Return calculated headloss over the pipe
  //! \retval headloss_ headloss over the pipe
  double headloss() const { return headloss_; }

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
