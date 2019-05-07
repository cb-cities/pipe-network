#include "node.h"

// compute polynomial coefficients for polynomial approximation of a function in
// a given interval a : third order,x^3, coefficient b : second_order,x^2,
// coefficient c : first order,x, coefficient d: zero order, C, coefficient.
Eigen::VectorXd pipenetwork::Node::compute_poly_coefficients(
    const std::array<double, 2>& x, const std::array<double, 2>& f,
    const std::array<double, 2>& df) {
  Eigen::VectorXd ret(4);
  double a =
      (2 * (f[0] - f[1]) - (x[0] - x[1]) * (df[1] + df[0])) /
      (std::pow(x[1], 3) - std::pow(x[0], 3) + 3 * x[0] * x[1] * (x[0] - x[1]));
  double b =
      (df[0] - df[1] + 3 * ((std::pow(x[1], 2) - std::pow(x[0], 2)) * a)) /
      (2 * (x[0] - x[1]));
  double c = df[1] - 3 * std::pow(x[1], 2) * a - 2 * x[1] * b;
  double d = f[1] - std::pow(x[1], 3) * a - std::pow(x[1], 2) * b - x[1] * c;
  ret << a, b, c, d;
  return ret;
}

Eigen::VectorXd pipenetwork::Node::get_pdd_poly_coef_1() {
  double x1 = minimum_pressure_;
  double f1 = 0;
  double x2 = minimum_pressure_ + pdd_smoothing_delta_;
  double f2 = std::pow(
      (x2 - minimum_pressure_) / (normal_pressure_ - minimum_pressure_), 0.5);
  double df1 = pdd_slope_;
  double df2 = 0.5 *
               std::pow((x2 - minimum_pressure_) /
                            (normal_pressure_ - minimum_pressure_),
                        -0.5) *
               (1 / (normal_pressure_ - minimum_pressure_));

  return compute_poly_coefficients({x1, x2}, {f1, f2}, {df1, df2});
}

Eigen::VectorXd pipenetwork::Node::get_pdd_poly_coef_2() {
  double x1 = normal_pressure_ - pdd_smoothing_delta_;
  double f1 = std::pow(
      (x1 - minimum_pressure_) / (normal_pressure_ - minimum_pressure_), 0.5);
  double x2 = normal_pressure_;
  double f2 = 1;
  double df1 = 0.5 *
               std::pow((x1 - minimum_pressure_) /
                            (normal_pressure_ - minimum_pressure_),
                        -0.5) *
               (1 / (normal_pressure_ - minimum_pressure_));
  double df2 = pdd_slope_;
  return compute_poly_coefficients({x1, x2}, {f1, f2}, {df1, df2});
}
