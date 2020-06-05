
#ifndef PIPE_NETWORK_PIPE_H
#define PIPE_NETWORK_PIPE_H
#include "link.h"

namespace pipenetwork {
//! Pipe Property
//! node1_name name of one end node
//! node2_name ame of the other end node
//! length, length of the pipe
//! diameter, diameter of the pipe
//! roughness, roughness of the pipe
//! status, status of the pipe (open or close)
struct PipeProp {
  std::string name{"None"};
  std::string node1_name{"None"};
  std::string node2_name{"None"};
  double length{std::numeric_limits<float>::max()};
  double diameter{std::numeric_limits<float>::max()};
  double roughness{std::numeric_limits<float>::max()};
  double minor_loss_coeff{0};
  LinkStatus status{LinkStatus::OPEN};
};

class Pipe : public Link {
 public:
  //! Constructor for a pipe
  //! \param[in] link_id link id
  //! \param[in] node1 one end node
  //! \param[in] node2 the other end node
  //! \param[in] pipe_prop struct with properties for the pipe
  Pipe(Index link_id, const Node& node1, const Node& node2,
       const PipeProp& pipe_prop)
      : Link(link_id, node1, node2), property_{pipe_prop} {};

  //! Virtual destructor
  ~Pipe() override{};

  //! Return pipe property
  PipeProp property() const { return property_; }

 private:
  // pipe properties
  PipeProp property_;
};

}  // namespace pipenetwork
#endif  // PIPE_NETWORK_PIPE_H
