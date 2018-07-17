#ifndef PIPE_NETWORK_NODE_H_
#define PIPE_NETWORK_NODE_H_

#include <Eigen/Dense>

//! Node class
//! \brief Class that stores the information about nodes
class Node {

 public:
  // Constructor with id and coordinates
  //! \param[in] id node id
  //! \param[in] coordinates coordinates of the node
  Node(unsigned id, const Eigen::Vector3d& coordinates)
      : id_{id}, coordinates_{coordinates} {}

  //! Destructor
  ~Node() = default;

  //! Copy constructor
  Node(const Node&) = delete;

  //! Assignment operator
  Node& operator=(const Node&) = delete;

  //! Move constructor
  Node(Node&&) = delete;

  //! Return id
  //! \retval id_ return id of the node
  unsigned id() const { return id_; }

  //! Return coordinates
  //! \retval coordinates_ return coordinates of the node
  Eigen::Vector3d coordinates() const { return coordinates_; }

  //! Return number of connection
  //! \retval nconnections_ return number of connection to the node
  unsigned nconnections() const { return nconnections_; }

  //! Assign hydraulic head at the node
  //! \param[in] head hydraulic head at the node
  void head(double head) {
    head_ = head;
    ishead_ = true;
  }

  //! Return hydraulic head
  //! \retval head_ return hydraulic head at the node
  double head() const { return head_; }

  //! Return head assignment status
  //! \retval ishead_ return head assignment status at the node
  bool ishead() const { return ishead_; }

  //! Assign discharge at the node
  //! \param[in] discharge discharge at the node
  void discharge(double discharge) {
    discharge_ = discharge;
    isdischarge_ = true;
  }

  //! Return discharge
  //! \retval discharge_ return discharge at the node
  double discharge() const { return discharge_; }

  //! Return discharge assignment status
  //! \retval isdischarge return discharge assignment status at the node
  bool isdischarge() const { return isdischarge_; }

 private:
  //! node id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! nodal coordinates
  Eigen::Vector3d coordinates_;
  //! number of connection to the node
  unsigned nconnections_{std::numeric_limits<unsigned>::max()};
  //! hydraulic head
  double head_{std::numeric_limits<double>::max()};
  //! discharge
  double discharge_{std::numeric_limits<double>::max()};
  //! whether head is assigned
  bool ishead_{false};
  //! whether discharge is assigned
  bool isdischarge_{false};
};

#endif  // PIPE_NETWORK_NODE_H_
