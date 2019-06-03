#ifndef PIPE_NETWORK_RESERVOIR_H
#define PIPE_NETWORK_RESERVOIR_H

#include "node_base.h"

namespace pipenetwork {

class Reservoir : public Node {
 public:
  //! Constructor with id abd head
  //! \param[in] id reservoir id
  //! \param[in] head base head for the reservoir
  Reservoir(Index id, const double head) : Node(id) {
    reservoir_info_["type"] = reservoir_type;
    reservoir_info_["head"] = head;
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
