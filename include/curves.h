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
  Pump_curve_prop(std::string curve_name,
                  std::vector<std::pair<double, double>>& curve_point);

  //! Name of the pump curve
  std::string name;
  //! Number of points used for curve generation
  unsigned num_points;
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
    // default pdd coefficients
    poly_coefficients["PDD_POLY_VEC1"] = {-18.749999999749996, 6.2499999999,
                                          1.000000082740371e-11,
                                          -4.440892098516782e-17};
    poly_coefficients["PDD_POLY_VEC2"] = {
        -0.6249920885505783, 37.249212823040864, -739.9780066609305,
        4900.811712406892};

    poly_coefficients["HW_POLY_VEC"] = {6619.952473405493, -2.562247355522429,
                                        0.0012305046454003125,
                                        3.4293453535907055e-09};
  };

  //! Destructor
  ~Curves() = default;
  //
  //! Method to add head pump curve
  //! \param[in] head_pump_curves a vector that contains information of each
  //! pump curve
  void add_pump_curves(const std::vector<Pump_curve_prop>& head_pump_props);

  std::map<std::string, Pump_curve_prop> pump_curves() const {
    return head_pump_curves;
  }
  std::map<std::string, Eigen::Vector4d> poly_coeffs() const {
    return poly_coefficients;
  }

 private:
  std::map<std::string, Eigen::Vector4d> poly_coefficients;
  std::map<std::string, Pump_curve_prop> head_pump_curves;
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
