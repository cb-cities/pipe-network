#ifndef PIPE_NETWORK_RESERVOIR_H
#define PIPE_NETWORK_RESERVOIR_H

#include "node.h"

namespace pipenetwork {
//! Reservoir Property
//! name junction name
//! head hydraulic head for the reservoir  for the junction
struct ReservoirProp {
  std::string name{"junction"};
  double head{std::numeric_limits<float>::min()};
};

class Reservoir : public Node {
 public:
  //! Constructor with id abd head
  //! \param[in] id reservoir id
  //! \param[in] head base head for the reservoir
  Reservoir(const Index id, const ReservoirProp res_prop)
      : Node(id), property_{res_prop} {};

  //! Return nodal info
  inline const ReservoirProp property() const { return property_; }

  //! Getter for head
  inline const double head() const { return property_.head; }

  //! Setter for head
  inline void head(const double head) { property_.head = head; }

  //! Discharge from the reservoir, mutable
  double discharge = 0;

 private:
  //! Reservoir information
  ReservoirProp property_;
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_RESERVOIR_H
