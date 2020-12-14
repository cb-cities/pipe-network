#ifndef PIPE_NETWORK_VALVE_H
#define PIPE_NETWORK_VALVE_H

namespace pipenetwork {
//! Valve Property
//! type valve type
//! setting pressure setting for PRV, PSV, or PBV, flow setting for FCV,  loss
//! coefficient for TCV,
struct ValveProp : LinkProp {
  ValveType type{ValveType::PRVALVE};
  double setting{0};
  double diameter{VALVE_DIAMETER};
  double minor_loss_coeff{MINOR_LOSS_COEFF};
};

//! Valve class
//! \brief Class that stores the information about pipes
class Valve : public Link {
 public:
  //! Constructor for a VALVE
  //! \param[in] link_id link id
  //! \param[in] node1 one end node
  //! \param[in] node2 the other end node
  //! \param[in] valve_prop struct with properties for the valve
  Valve(Index link_id, const Node& node1, const Node& node2,
        const ValveProp& valve_prop)
      : Link(link_id, node1, node2), property_{valve_prop} {};

  //! Virtual destructor
  ~Valve() override{};

  //! Return Valve property
  const ValveProp& property() const { return property_; }

 private:
  //! valve properties
  ValveProp property_;
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_VALVE_H
