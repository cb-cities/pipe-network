#ifndef PIPE_NETWORK_PUMP_H
#define PIPE_NETWORK_PUMP_H

namespace pipenetwork {

//! Pump Property
//! node1 ptr for end node
//! node2 ptr for end node
//! type pump type
//! curve_name pump head curve name
//! speed speed for the pump
//! pattern pattern for speed setting
struct Pump_prop {
  std::string id;
  std::shared_ptr<pipenetwork::Node> node1{NULL};
  std::shared_ptr<pipenetwork::Node> node2{NULL};
  Link_type pump_type{POWERPUMP};
  Link_status pump_status{OPEN};
  int curve_name{-1};
  double power{50};
  double speed{1.0};
  std::string node1_id{"None"};
  std::string node2_id{"None"};
  std::string Pattern{"None"};
};

class Pump : public Link {
 public:
  //! Constructor with pump property, which contains all the information for
  //! the valve \param[in] pipe_prop struct with properties for the pipe
  Pump(const Pump_prop& pump_prop)
      : Link(pump_prop.id, pump_prop.node1, pump_prop.node2,pump_prop.pump_status) {
    pump_info_["type"] = pump_prop.pump_type;
    pump_info_["curve_name"] = pump_prop.curve_name;
    pump_info_["power"] = pump_prop.power;
    pump_info_["speed"] = pump_prop.speed;
  };

  //! Virtual destructor
  ~Pump() override{};

  //! Return link info
  std::map<std::string, double> link_info() const override {
    return pump_info_;
  }

 private:
  //! valve information
  std::map<std::string, double> pump_info_;
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_PUMP_H
