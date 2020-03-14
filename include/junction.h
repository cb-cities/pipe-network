#ifndef PIPE_NETWORK_JUNCTION_H
#define PIPE_NETWORK_JUNCTION_H

#include "node_base.h"

namespace pipenetwork {
//! Junction Property
//! id junction id
//! elevation elevation for the junction
//! demand base demand for the junction
//! leak_diameter diameter of the leak hole for the junction
struct Junction_prop {
  std::string id;
  double elevation{std::numeric_limits<float>::min()};
  double demand{std::numeric_limits<float>::min()};
  double leak_diameter{0};
};

class Junction : public Node {
 public:
  //! Constructor with id, elevation, demand and leak diameter
  //! \param[in] junc_prop junction properties

  Junction(const Junction_prop& junc_prop) : Node(junc_prop.id) {

    junction_info_["type"] = JUNCTION;
    junction_info_["elevation"] = junc_prop.elevation;
    junction_info_["demand"] = junc_prop.demand;
    junction_info_["leak_area"] =
        std::pow((junc_prop.leak_diameter / 2), 2) * PI;

    update_sim_demand(junc_prop.demand);
    update_sim_head(junc_prop.elevation);
  };

  //! Virtual destructor
  ~Junction() override{};

  //! Delete Copy constructor
  Junction(const Junction&) = delete;

  //! Delete Assignment operator
  Junction& operator=(const Junction&) = delete;

  //! Move constructor
  Junction(Junction&&) = delete;

  //! Return nodal info
  std::map<std::string, double> nodal_info() const override {
    return junction_info_;
  }

 private:
  // node information, has key : type, elevation, demand, leak_area
  std::map<std::string, double> junction_info_;
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_JUNCTION_H
