#ifndef PIPE_NETWORK_PUMP_H
#define PIPE_NETWORK_PUMP_H

#include "link.h"

namespace pipenetwork {
//! Pump Property
//! type pump type
//! curve_name pump head curve id
//! speed speed for the pump
//! pattern pattern for speed setting
struct PumpProp : LinkProp {
  PumpType type{PumpType::POWERPUMP};
  int curve_id{-1};
  double power{PUMP_POWER};
  double speed{PUMP_SPEED};
};

//! Pump class
//! \brief Class that stores the information about pumps
class Pump : public Link {
 public:
  //! Constructor with pump property, which contains all the information for
  //! the valve \param[in] pipe_prop struct with properties for the pipe
  Pump(Index link_id, const Node& node1, const Node& node2,
       const PumpProp& pump_prop)
      : Link(link_id, node1, node2), property_{pump_prop} {};

  //! Virtual destructor
  ~Pump() override{};

  //! Return pump property
  const PumpProp& property() const { return property_; }

 private:
  //! valve information
  PumpProp property_;
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_PUMP_H
