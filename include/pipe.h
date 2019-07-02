
#ifndef PIPE_NETWORK_PIPE_H
#define PIPE_NETWORK_PIPE_H
#include "link.h"

namespace pipenetwork {
//! Pipe Property
//! node1 ptr for end node
//! node2 ptr for end node
//! length, length of the pipe
//! diameter, diameter of the pipe
//! roughness, roughness of the pipe
//! status, status of the pipe (open or close)
struct Pipe_prop {
  std::string id;
  std::shared_ptr<pipenetwork::Node> node1{NULL};
  std::shared_ptr<pipenetwork::Node> node2{NULL};
  double length{std::numeric_limits<float>::max()};
  double diameter{std::numeric_limits<float>::max()};
  double roughness{std::numeric_limits<float>::max()};
  Pipe_status status{OPEN};
  std::string node1_id{"None"};
  std::string node2_id{"None"};
};

class Pipe : public Link {
 public:
  //! Constructor with two end nodes, length, diameter roughness and pipe status
  //! \param[in] pipe_prop struct with properties for the pipe
  Pipe(const Pipe_prop& pipe_prop)
      : Link(pipe_prop.id, pipe_prop.node1, pipe_prop.node2) {
    pipe_info_["type"] = PIPE;
    pipe_info_["length"] = pipe_prop.length;
    pipe_info_["diameter"] = pipe_prop.diameter;
    pipe_info_["roughness"] = pipe_prop.roughness;
    pipe_info_["status"] = pipe_prop.status;
  };

  //! Virtual destructor
  ~Pipe() override{};

  //! Return link info
  std::map<std::string, double> link_info() const override {
    return pipe_info_;
  }

 private:
  // node information, has key : type, elevation, demand, leak_area
  std::map<std::string, double> pipe_info_;
};

}  // namespace pipenetwork
#endif  // PIPE_NETWORK_PIPE_H
