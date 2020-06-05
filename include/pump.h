#ifndef PIPE_NETWORK_PUMP_H
#define PIPE_NETWORK_PUMP_H

#include "link.h"

namespace pipenetwork {
//! Pump Property
//! node1_name name of one end node
//! node2_name ame of the other end node
//! type pump type
//! curve_name pump head curve id
//! speed speed for the pump
//! pattern pattern for speed setting
struct PumpProp {
  std::string name{"None"};
  std::string node1_name{"None"};
  std::string node2_name{"None"};
  PumpType type{PumpType::POWERPUMP};
  LinkStatus status{LinkStatus::OPEN};
  int curve_id{-1};
  double power{PUMP_POWER};
  double speed{PUMP_SPEED};
};

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
  PumpProp property() const { return property_; }

 private:
  //! valve information
  PumpProp property_;
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_PUMP_H
