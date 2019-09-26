#ifndef PIPE_NETWORK_CURVES_H
#define PIPE_NETWORK_CURVES_H

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "settings.h"

namespace pipenetwork {
//! Pump_curve_prop struct
//! \brief Struct for pump curves that contains curve information for a given
//! head pump curve
struct Pump_curve_prop {
  Pump_curve_prop() = default;
  //! Constructor
  //! \param[in] curve_name name of the curve
  //! \param[in] curve_point vector (1 or 3) of curve points (flow, head) used
  //! for pump curve construction
  Pump_curve_prop(std::string& curve_name,
                  std::vector<std::pair<double, double>>& curve_point);

  //! Name of the pump curve
  std::string name;
  //! Number of points used for curve generation
  unsigned num_points{0};
  //! Points used for curve generation
  std::vector<std::pair<double, double>> points;
  //! Head curve coefficients: H = head_curve_coefficients[0] -
  //! head_curve_coefficients[1]*Q^head_curve_coefficients[2]
  Eigen::Vector3d head_curve_coefficients{0, 0, 0};
  //! Pump polynomial approximation coefficients
  Eigen::Vector4d poly_coefficients{0, 0, 0, 0};
  //! Pump line approximation coefficients
  Eigen::Vector2d line_param{0, 0};

  //! Method to compute pump polynomial approximation coefficients
  Eigen::Vector4d get_pump_poly_coefficients();
  //! Method to compute pump line approximation coefficients
  Eigen::Vector2d get_pump_line_params();
};

//! Curves Class
//! \brief Class to hold all the curves information in hydraulic simulation
class Curves {
 public:
  //! Constructor
  Curves() {
    // default pdd(pressure demand driven) and hw(harzian williams) coefficients
    poly_coefficients_["PDD_POLY_VEC1"] = compute_pdd1_poly_vec();
    poly_coefficients_["PDD_POLY_VEC2"] = compute_pdd2_poly_vec();
    poly_coefficients_["HW_POLY_VEC"] = compute_hw_poly_vec();
  };

  //! Destructor
  ~Curves() = default;

  //
  //! Method to add head pump curve
  //! \param[in] head_pump_curves a vector that contains information of each
  //! pump curve
  void add_pump_curves(const std::vector<Pump_curve_prop>& head_pump_props);

  //! Method to add polynomial coefficients for leak equation
  //! \param[in] node_name name of the leaky node
  //! \param[in] leak_area area of the leak hole
  void add_leak_poly_vec(std::string& node_name, double leak_area) {
    poly_coefficients_[node_name] = compute_leak_poly_vec(leak_area);
  };

  //! get the head pump curves map
  std::map<std::string, Pump_curve_prop> pump_curves() const {
    return head_pump_curves_;
  }

  //! get the poly coefficient curves map
  std::map<std::string, Eigen::Vector4d> poly_coeffs() const {
    return poly_coefficients_;
  }

  //! get pump string name to int code map
  int pump_str_int(std::string k) const { return pump_str_int_.at(k); }

  //! get pump int code to string name map
  std::string pump_int_str(int k) const { return pump_int_str_.at(k); }

 private:
  //! map that stores all the polynomial approximation coefficients
  std::map<std::string, Eigen::Vector4d> poly_coefficients_;
  //! map that stores all the head pump curves information
  std::map<std::string, Pump_curve_prop> head_pump_curves_;
  //! map that convert pump string key to int
  std::map<std::string, int> pump_str_int_;
  std::map<int, std::string> pump_int_str_;
  //! Method to compute the polynomial coefficients for harzian williams
  //! equation
  //! \retval poly_coef polynomial approximation coefficient
  Eigen::Vector4d compute_hw_poly_vec();
  //! Method to compute the polynomial coefficients for the first part of
  //! pressure demand driven demand equation
  //! \retval poly_coef polynomial approximation coefficient
  Eigen::Vector4d compute_pdd1_poly_vec();
  //! Method to compute the polynomial coefficients for the second part of
  //! pressure demand driven demand equation
  //! \retval poly_coef polynomial approximation coefficient
  Eigen::Vector4d compute_pdd2_poly_vec();
  //! Method to compute the polynomial coefficients for leak equations
  //! \param[in] leak_area leak area of the leak hole
  //! \retval leak_poly_coef polynomial approximation coefficient for leak
  //! equation
  Eigen::Vector4d compute_leak_poly_vec(double leak_area);
};

//! Method to compute the coefficients of a smoothing polynomial
//! \param[in] x points on the x-axis at which the smoothing polynomial begins
//! and ends
//! \param[in] f function evaluated at x1 and f2 \param[in] df
//! derivative evaluated at x1 and x2
//! \retval  A vector with the smoothing
//! polynomail coefficients starting with the cubic term.
Eigen::Vector4d compute_poly_coefficients(const std::array<double, 2>& x,
                                          const std::array<double, 2>& f,
                                          const std::array<double, 2>& df);
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_CURVES_H
