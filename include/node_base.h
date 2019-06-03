#ifndef PIPE_NETWORK_NODE_H_
#define PIPE_NETWORK_NODE_H_

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <map>
#include <string>
#include <memory>

#include "settings.h"

namespace pipenetwork {

//! Node class
//! \brief Base Class that stores the information about nodes
class Node {

 public:
  Node(Index node_id) : node_id_{node_id} {};

  //! Destructor
  virtual ~Node(){};

  //! Delete Copy constructor
  Node(const Node&) = delete;

  //! Delete Assignment operator
  Node& operator=(const Node&) = delete;

  //! Move constructor
  Node(Node&&) = delete;

  //! Return nodal info
  virtual std::map<std::string, double> nodal_info() const = 0;

  //! Return nodal id
  Index id() const { return node_id_; }

  //! Assign simulated demand
  void sim_demand(double demand) { sim_demand_ = demand; }

  //! Return simulated demand
  double sim_demand() const { return sim_demand_; }

  //! Assign simulated head
  void sim_head(double head) { sim_head_ = head; }

  //! Return simulated head
  double sim_head() const { return sim_head_; }

  //! Assign simulated leak discharge
  void sim_leak(double leak) { sim_leak_ = leak; }

  //! Return simulated leak discharge
  double sim_leak() const { return sim_leak_; }

 private:
  Index node_id_;
  // demand from simulation
  double sim_demand_{0};
  // head from simulation
  double sim_head_{0};
  // leak discharge from simulation
  double sim_leak_{0};
};
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_NODE_H_
