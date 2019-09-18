#ifndef PIPE_NETWORK_VALVE_H
#define PIPE_NETWORK_VALVE_H

namespace pipenetwork {
//! Valve Property
//! node1 ptr for end node
//! node2 ptr for end node
//! type valve type
//! setting pressure setting for PRV, PSV, or PBV, flow setting for FCV,  loss
//! coefficient for TCV,
struct Valve_prop {
  std::string id;
  std::shared_ptr<pipenetwork::Node> node1{NULL};
  std::shared_ptr<pipenetwork::Node> node2{NULL};
  Link_type valve_type{PRVALVE};
  Link_status valve_status{ACTIVE};
  double setting{0};
  double diameter{0.3048};
  double minor_loss{0};
  std::string node1_id{"None"};
  std::string node2_id{"None"};
};

class Valve : public Link {
 public:
  //! Constructor with valve property, which contains all the information for
  //! the valve
  //! \param[in] pipe_prop struct with properties for the pipe
  explicit Valve(const Valve_prop& valve_prop)
      : Link(valve_prop.id, valve_prop.node1, valve_prop.node2,valve_prop.valve_status) {
    valve_info_["type"] = valve_prop.valve_type;
    valve_info_["setting"] = valve_prop.setting;
    valve_info_["minor_loss"] = valve_prop.minor_loss;
    valve_info_["diameter"] = valve_prop.diameter;
  };

  //! Virtual destructor
  ~Valve() override{};

  //! Return link info
  std::map<std::string, double> link_info() const override {
    return valve_info_;
  }

 private:
  //! valve information
  std::map<std::string, double> valve_info_;
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_VALVE_H
