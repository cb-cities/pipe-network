#ifndef PIPE_NETWORK_JUNCTION_H
#define PIPE_NETWORK_JUNCTION_H

#include "node.h"
namespace pipenetwork {
//! Junction Property
//! name junction name
//! elevation elevation for the junction
//! demand base demand for the junction
//! leak_diameter diameter of the leak hole for the junction
struct JunctionProp {
  std::string name{"junction"};
  double elevation{std::numeric_limits<float>::min()};
  double demand{std::numeric_limits<float>::min()};
  double leak_diameter{0};
};

class Junction : public Node {
 public:
  //! Constructor with id, elevation, demand and leak diameter
  //! \param[in] junc_prop junction properties

  Junction(const Index id, const JunctionProp& junc_prop)
      : Node(id), property_{junc_prop} {};

  //! Return nodal info
  inline const JunctionProp property() const { return property_; }

  //! Return leak area
  inline const double leak_area() const {
    return std::pow((property_.leak_diameter / 2), 2) * PI;
  }

  //! Junction demand, mutable
  double demand{0};
  //! Junction head, mutable
  double head{1e-3};
  //! Junction leak discharge, mutable
  double leak_discharge{0};

 private:
  //! junction information
  JunctionProp property_;
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_JUNCTION_H
