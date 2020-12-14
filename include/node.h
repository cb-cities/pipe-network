#ifndef PIPE_NETWORK_NODE_H_
#define PIPE_NETWORK_NODE_H_

#include <string>
#include <variant>

#include "settings.h"

namespace pipenetwork {

//! Node class
//! \brief Base Class that stores the information about nodes
class Node {
 public:
  explicit Node(const Index node_id) : node_id_{node_id} {};

  //! Destructor
  virtual ~Node(){};

  //! Return nodal id
  Index id() const { return node_id_; }

 private:
  Index node_id_;
};
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_NODE_H_
