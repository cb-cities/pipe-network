
#ifndef PIPE_NETWORK_PIPE_H
#define PIPE_NETWORK_PIPE_H
#include "link.h"

namespace pipenetwork {

class Pipe : public Link {
 public:
  //! Constructor with two end nodes, length, diameter roughness and pipe status
  //! \param[in] node1 ptr for end node
  //! \param[in] node2 ptr for end node
  //! \param[in] length, length of the pipe
  //! \param[in] diameter, diameter of the pipe
  //! \param[in] roughness, roughness of the pipe
  //! \param[in] status, status of the pipe (open or close)
  Pipe(Index id, const std::shared_ptr<pipenetwork::Node>& node1,
       const std::shared_ptr<pipenetwork::Node>& node2, const double length,
       const double diameter, const double roughness, const Pipe_status status)
      : Link(id, node1, node2) {
    pipe_info_["type"] = PIPE;
    pipe_info_["length"] = length;
    pipe_info_["diameter"] = diameter;
    pipe_info_["roughness"] = roughness;
    pipe_info_["status"] = status;
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
