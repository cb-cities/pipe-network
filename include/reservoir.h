#ifndef PIPE_NETWORK_RESERVOIR_H
#define PIPE_NETWORK_RESERVOIR_H

#include "node_base.h"

//! Reservoir Property
//! id junction id
//! head head for the reservoir
namespace pipenetwork {
struct Reservoir_prop {
  std::string id;
  double head{std::numeric_limits<float>::min()};
};

class Reservoir : public Node {
 public:
  //! Constructor with id abd head
  //! \param[in] id reservoir id
  //! \param[in] head base head for the reservoir
  Reservoir(const Reservoir_prop& res_prop) : Node(res_prop.id) {
    reservoir_info_["type"] = RESERVOIR;
    reservoir_info_["head"] = res_prop.head;
    update_sim_demand(0);
    update_sim_head(res_prop.head);
  };

  //! Virtual destructor
  ~Reservoir() override{};

  //! Delete Copy constructor
  Reservoir(const Reservoir&) = delete;

  //! Delete Assignment operator
  Reservoir& operator=(const Reservoir&) = delete;

  //! Move constructor
  Reservoir(Reservoir&&) = delete;

  //! Return nodal info
  std::map<std::string, double> nodal_info() const override {
    return reservoir_info_;
  }

 private:
  //! node information, has key : type, elevation, demand, leak_area
  std::map<std::string, double> reservoir_info_;
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_RESERVOIR_H
