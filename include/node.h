#ifndef PIPE_NETWORK_NODE_H_
#define PIPE_NETWORK_NODE_H_

#include <Eigen/Dense>

namespace pipenetwork {

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
  //! \retval id_ id of the node
  unsigned id() const { return id_; }

  //! Return coordinates
  //! \retval coordinates_ coordinates of the node
  Eigen::Vector3d coordinates() const { return coordinates_; }

  //! Return number of connection
  //! \retval nconnections_ number of connection to the node
  unsigned nconnections() const { return nconnections_; }

  //! Assign hydraulic head at the node
  //! \param[in] head hydraulic head at the node
  void head(double head) {
    head_ = head;
    iter_head_ = head;
    ishead_ = true;
  }
  void iter_head(double iter_head) { iter_head_ = iter_head; }

  bool isres() const { return isres_; }

  //! Return hydraulic head
  //! \retval head_ hydraulic head at the node
  double head() const { return head_; }
  double iter_head() const { return iter_head_; }

  //! Return head assignment status
  //! \retval ishead_ head assignment status at the node
  bool ishead() const { return ishead_; }

  //! Assign discharge at the node
  //! \param[in] discharge discharge at the node
  void discharge(double discharge) {
    discharge_ = discharge;
    iter_discharge_ = discharge;
    isdischarge_ = true;
    // if demand is negative, it is a reservoir or tank
    if (discharge < 0) {
      isres_ = true;
    }
  }

  double iter_discharge() { return iter_discharge_; }

  void iter_discharge(double iter_discharge) {
    iter_discharge_ = iter_discharge;
  }
  //! Return discharge
  //! \retval discharge_ discharge at the node
  double discharge() const { return discharge_; }

  //! Return discharge assignment status
  //! \retval isdischarge_ discharge assignment status at the node
  bool isdischarge() const { return isdischarge_; }

  //! min pressure
  void min_pressure(double min_pressure){
      minimum_pressure_ = min_pressure;
  }
  double min_pressure() const { return  minimum_pressure_;}

  //! norm pressure
  void norm_pressure(double norm_pressure){
      normal_pressure_ = norm_pressure;
  }
  double norm_pressure() const {return normal_pressure_;};

  //!pdd delta
  double pdd_smooth_delta() const {return pdd_smoothing_delta_;}

  //!pdd slope
  double pdd_slope() const {return pdd_slope_;}



  Eigen::VectorXd get_pdd_poly_coef_1();
    Eigen::VectorXd get_pdd_poly_coef_2();



 private:
  //! node id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! nodal coordinates
  Eigen::Vector3d coordinates_;
  //! number of connection to the node
  unsigned nconnections_{std::numeric_limits<unsigned>::max()};
  //! hydraulic head
  double head_{0};
  double iter_head_{50};
  //! discharge
  double discharge_{10};
  double iter_discharge_{9};
  //! whether head is assigned
  bool ishead_{false};
  //! whether discharge is assigned
  bool isdischarge_{false};
  //! whether node is reservoir
  bool isres_{false};
  //! minimum_pressure
  double minimum_pressure_{0};
  //! normal pressure
  double normal_pressure_{20};
  //! pdd smoothing delta
  double pdd_smoothing_delta_{0.2};
  double pdd_slope_{1e-12};

    //! Method to compute the coefficients of a smoothing polynomial
    //! \param[in] x1 point on the x-axis at which the smoothing polynomial begins
    //! \param[in] x2 point on the x-axis at which the smoothing polynomial ends
    //! \param[in] f1 function evaluated at x1
    //! \param[in] f2 function evaluated at x2
    //! \param[in] df1 derivative evaluated at x1
    //! \param[in] df2 derivative evaluated at x2
    //! \retval  A vector with the smoothing polynomail coefficients starting with
    //! the cubic term.
    Eigen::VectorXd compute_poly_coefficients(double x1, double x2, double f1,
                                              double f2, double df1, double df2);
};
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_NODE_H_
