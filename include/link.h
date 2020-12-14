#ifndef PIPE_NETWORK_LINK_H
#define PIPE_NETWORK_LINK_H

#include <array>
#include <memory>
#include <string>

#include "node.h"

namespace pipenetwork {
//! Link Property
//! name name of the link
//! node1_name name of one end node
//! node2_name ame of the other end node
struct LinkProp {
  std::string name{"None"};
  std::string node1_name{"None"};
  std::string node2_name{"None"};
  LinkStatus status{LinkStatus::OPEN};
};

//! Link class
//! \brief Base Class that stores the information about links
class Link {
 public:
  //! Constructor
  //! \param[in] link_id link id
  //! \param[in] node1 one end node
  //! \param[in] node2 the other end node
  Link(const Index link_id, const Node& node1, const Node& node2)
      : id_{link_id} {
    nodes_ = std::make_pair(std::make_shared<Node>(node1),
                            std::make_shared<Node>(node2));
  };

  //! Destructor
  virtual ~Link(){};

  //! Return link id
  Index id() const { return id_; }

  //! Return end nodes
  const std::pair<std::shared_ptr<Node>, std::shared_ptr<Node>>& nodes() const {
    return nodes_;
  }

 private:
  //! link id
  Index id_;
  //! pair of node pointers which form the link
  std::pair<std::shared_ptr<Node>, std::shared_ptr<Node>> nodes_;
};
}  // namespace pipenetwork
;

#endif  // PIPE_NETWORK_LINK_H
