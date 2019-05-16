#ifndef PIPE_NETWORK_NODE_H_
#define PIPE_NETWORK_NODE_H_

#include <Eigen/Dense>
#include <array>
#include <cmath>

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
  void elevation(double head) {
    elevation_ = head;
    head_ = head;
    ishead_ = true;
  }
  void head(double iter_head) { head_ = iter_head; }

  //! Return hydraulic head
  //! \retval head_ hydraulic head at the node
  double elevation() const { return elevation_; }
  double head() const { return head_; }

  //! Return if node is a reservoir
  //! \retval isres_ reservoir status
  bool isres() const { return isres_; }

  //! Assign discharge at the node
  //! \param[in] discharge discharge at the node
  void demand(double discharge) {
    demand_ = discharge;
    iter_demand_ = discharge;
    isdischarge_ = true;
    // if demand is negative, it is a reservoir or tank
    if (discharge < 0) {
      isres_ = true;
    }
  }
  double iter_demand() { return iter_demand_; }

  void iter_demand(double iter_discharge) { iter_demand_ = iter_discharge; }
  //! Return discharge
  //! \retval discharge_ discharge at the node
  double demand() const { return demand_; }

  //! Return discharge assignment status
  //! \retval isdischarge_ discharge assignment status at the node
  bool isdischarge() const { return isdischarge_; }

  //! min pressure
  double min_pressure() const { return minimum_pressure_; }

  //! norm pressure
  double norm_pressure() const { return normal_pressure_; };

  //! pdd delta
  double pdd_smooth_delta() const { return pdd_smoothing_delta_; }
  //! pdd slope
  double pdd_slope() const { return pdd_slope_; }

  //! leak status
  bool is_leak() const {
      return leak_area_ > 0 || false;
  }
  //! set diameter/area for the leak
  void leak_diameter(double leak_diameter) {
    leak_diameter_ = leak_diameter;
    leak_area_ = std::pow((leak_diameter_ / 2), 2) * pi_;
  }

  double leak_area() const { return leak_area_; }
  double leak_discharge() const { return leak_diameter_; }
  double iter_leak_discharge() { return iter_leak_discharge_; }

  void iter_leak_discharge(double iter_leak_discharge) {
    iter_leak_discharge_ = iter_leak_discharge;
  }

  double g() const { return g_; }
  double leak_discharge_coefficient() const {
    return leak_discharge_coefficient_;
  }

  //! polynomial coefficients for pressure demand equations
  Eigen::VectorXd get_pdd_poly_coef_1();
  Eigen::VectorXd get_pdd_poly_coef_2();
  //! polynomial coefficients for leakage equations
  Eigen::VectorXd get_leak_poly_coef();

 private:
  //! node id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! nodal coordinates
  Eigen::Vector3d coordinates_;
  //! number of connection to the node
  unsigned nconnections_{std::numeric_limits<unsigned>::max()};
  //! hydraulic head
  double elevation_{0};
  double head_{0};
  //! discharge
  double demand_{10};
  double iter_demand_{10};
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

  //! leak diameter
  double leak_diameter_{0};
  //! Value of pi
  double pi_{3.1415926};
  //! leak area
  double leak_area_{0};

  //! accleration due to gravity
  double g_{9.81};
  //! leak dischage (static and iter)
  double leak_discharge_{0};
  double iter_leak_discharge_{0};
  //! leak coefficients
  double leak_discharge_coefficient_{0.75};

  //! polynomial coefficients for demand-pressure leak-pressure
  bool is_pdd_coef1_{false};
  bool is_pdd_coef2_{false};
  bool is_leak_coef_{false};
  Eigen::VectorXd pdd_poly_coef_1_;
  Eigen::VectorXd pdd_poly_coef_2_;
  Eigen::VectorXd leak_poly_coef_;

  //! Method to compute the coefficients of a smoothing polynomial
  //! \param[in] x points on the x-axis at which the smoothing polynomial begins
  //! and ends \param[in] f function evaluated at x1 and f2 \param[in] df
  //! derivative evaluated at x1 and x2 \retval  A vector with the smoothing
  //! polynomail coefficients starting with the cubic term.
  Eigen::VectorXd compute_poly_coefficients(const std::array<double, 2>& x,
                                            const std::array<double, 2>& f,
                                            const std::array<double, 2>& df);
  //! Method to compute the polynomial coefficients for pressure demand
  //! equations
  void compute_pdd_poly_coef_1();
  void compute_pdd_poly_coef_2();
  //! Method to compute the polynomial coefficients for leak equations
  void compute_leak_poly_coef();
};
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_NODE_H_
