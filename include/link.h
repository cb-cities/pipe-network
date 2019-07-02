#ifndef PIPE_NETWORK_LINK_H
#define PIPE_NETWORK_LINK_H

#include <Eigen/Dense>
#include <array>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>

#include "node_base.h"
#include "settings.h"

namespace pipenetwork {

//! Link class
//! \brief Base Class that stores the information about links
class Link {

 public:
  Link(std::string link_id, const std::shared_ptr<pipenetwork::Node>& node1,
       const std::shared_ptr<pipenetwork::Node>& node2)
      : link_id_{link_id} {
    if (!node1 || !node2)
      throw std::invalid_argument(
          "Link can not be constructed due to invalid input node pointers");
    nodes_ = std::make_pair(node1, node2);
  };

  //! Destructor
  virtual ~Link(){};

  //! Delete Copy constructor
  Link(const Link&) = delete;

  //! Delete Assignment operator
  Link& operator=(const Link&) = delete;

  //! Move constructor
  Link(Link&&) = delete;

  //! Return link info
  virtual std::map<std::string, double> link_info() const = 0;

  //! Return link id
  std::string id() const { return link_id_; }
  //! Return end nodes
  std::pair<std::shared_ptr<pipenetwork::Node>,
            std::shared_ptr<pipenetwork::Node>>
      nodes() const {
    return nodes_;
  }

  //! Assign simulated discharge
  void update_sim_discharge(double sim_demand) { sim_demand_ = sim_demand; }

  //! Return simulated discharge
  double sim_discharge() const { return sim_demand_; }

 private:
  //! link id
  std::string link_id_;
  //! pair of node pointers which form the pipe
  std::pair<std::shared_ptr<pipenetwork::Node>,
            std::shared_ptr<pipenetwork::Node>>
      nodes_;
  //! demand from simulation
  double sim_demand_{0};
};
}  // namespace pipenetwork
;

#endif  // PIPE_NETWORK_LINK_H
