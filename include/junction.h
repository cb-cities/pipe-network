#ifndef PIPE_NETWORK_JUNCTION_H
#define PIPE_NETWORK_JUNCTION_H

#include "node_base.h"

namespace pipenetwork {

class Junction : public Node {
 public:
  //! Constructor with id, elevation, demand and leak diameter
  //! \param[in] id junction id
  //! \param[in] elevation elevation for the junction
  //! \param[in] demand base demand for the junction
  //! \param[in] leak_diameter diameter of the leak hole for the junction
  Junction(Index id, const double elevation, const double demand,
           const double leak_diameter)
      : Node(id) {
    junction_info_["type"] = junction_type;
    junction_info_["elevation"] = elevation;
    junction_info_["demand"] = demand;
    junction_info_["leak_area"] = std::pow((leak_diameter / 2), 2) * PI;
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
