
#include "curves.h"

// Returns the A, B, C coefficients for a 1-point or a 3-point pump curve.
// Coefficient can only be calculated for pump curves.
//
// For a single point curve the coefficients are generated according to the
// following equation:
//
// A = 4/3 * H_1
// B = 1/3 * H_1/Q_1^2
// C = 2
//
// For a three point curve the coefficients are generated according to the
// following equation: When the first point is a zero flow: (All INP files we
// have come across)
//
// A = H_1
// C = ln((H_1 - H_2)/(H_1 - H_3))/ln(Q_2/Q_3)
// B = (H_1 - H_2)/Q_2^C
//
// When the first point is not zero, numpy fsolve is called to solve the
// following system of equation:
//
// H_1 = A - B*Q_1^C
// H_2 = A - B*Q_2^C
// H_3 = A - B*Q_3^C
//
// Multi point curves are currently not supported
pipenetwork::Pump_curve_prop::Pump_curve_prop(
    std::string& curve_name, std::vector<std::pair<double, double>> curve_point)
    : name{curve_name}, points{curve_point} {
  num_points = points.size();
  double A, B, C;
  if (num_points == 1) {
    auto Q_1 = points[0].first;
    auto H_1 = points[0].second;

    A = (4.0 / 3.0) * H_1;
    B = (1.0 / 3.0) * (H_1 / (std::pow(Q_1, 2)));
    C = 2;
  } else if (num_points == 3) {
    auto Q_1 = points[0].first;
    auto H_1 = points[0].second;
    auto Q_2 = points[1].first;
    auto H_2 = points[1].second;
    auto Q_3 = points[2].first;
    auto H_3 = points[2].second;

    if (Q_1 == 0.0) {
      A = H_1;
      C = std::log((H_1 - H_2) / (H_1 - H_3)) / std::log(Q_2 / Q_3);
      B = (H_1 - H_2) / (std::pow(Q_2, C));
    }
    //! TODO: HEAD_CURVE_COEFFICIENTS FOR 3 POINTS WITH NONZERO STARTING POINT
  } else {
    throw std::runtime_error(
        "Coefficient for Multipoint pump curves cannot be generated");
  }
  if (A < 0 || B < 0 || C < 0) {
    throw std::runtime_error(
        "Value of pump head curve coefficient is negative, which is not "
        "allowed");
  }
  head_curve_coefficients = {A, B, C};
  if (C < 1) {
    poly_coefficients = get_pump_poly_coefficients();
  } else {
    line_param = get_pump_line_params();
  }
}

Eigen::Vector4d pipenetwork::Pump_curve_prop::get_pump_poly_coefficients() {

  auto f1 = PUMP_M * PUMP_Q1 + head_curve_coefficients[0];
  auto f2 = head_curve_coefficients[0] -
            head_curve_coefficients[1] *
                std::pow(PUMP_Q2, head_curve_coefficients[2]);
  auto df1 = PUMP_M;
  auto df2 = -head_curve_coefficients[1] * head_curve_coefficients[2] *
             std::pow(PUMP_Q2, (head_curve_coefficients[2] - 1.0));

  return compute_poly_coefficients({PUMP_Q1, PUMP_Q2}, {f1, f2}, {df1, df2});
}

Eigen::Vector2d pipenetwork::Pump_curve_prop::get_pump_line_params() {
  auto q_bar = std::pow(
      (PUMP_M / (-head_curve_coefficients[1] * head_curve_coefficients[2])),
      (1.0 / (head_curve_coefficients[2] - 1.0)));
  auto h_bar =
      head_curve_coefficients[0] -
      head_curve_coefficients[1] * std::pow(q_bar, head_curve_coefficients[2]);
  return {q_bar, h_bar};
}

Eigen::Vector4d pipenetwork::compute_poly_coefficients(
    const std::array<double, 2>& x, const std::array<double, 2>& f,
    const std::array<double, 2>& df) {
  Eigen::Vector4d ret;
  double a =
      (2 * (f[0] - f[1]) - (x[0] - x[1]) * (df[1] + df[0])) /
      (std::pow(x[1], 3) - std::pow(x[0], 3) + 3 * x[0] * x[1] * (x[0] - x[1]));
  double b =
      (df[0] - df[1] + 3 * ((std::pow(x[1], 2) - std::pow(x[0], 2)) * a)) /
      (2 * (x[0] - x[1]));
  double c = df[1] - 3 * std::pow(x[1], 2) * a - 2 * x[1] * b;
  double d = f[1] - std::pow(x[1], 3) * a - std::pow(x[1], 2) * b - x[1] * c;
  ret = {a, b, c, d};
  return ret;
}

void pipenetwork::Curves::add_pump_curves(
    const std::vector<Pump_curve_prop>& head_pump_props) {
  int count = 0;
  for (const auto& pump_prop : head_pump_props) {
    head_pump_curves_[pump_prop.name] = pump_prop;
    pump_str_int_[pump_prop.name] = count;
    pump_int_str_[count] = pump_prop.name;
    ++count;
  }
}

Eigen::Vector4d pipenetwork::Curves::compute_hw_poly_vec() {
  auto x1 = HW_Q1;
  auto x2 = HW_Q2;
  auto f1 = HW_M * HW_Q1;
  auto f2 = std::pow(HW_Q2, 1.852);
  auto df1 = HW_M;
  auto df2 = 1.852 * std::pow(HW_Q2, 0.852);
  return compute_poly_coefficients({x1, x2}, {f1, f2}, {df1, df2});
}

Eigen::Vector4d pipenetwork::Curves::compute_pdd1_poly_vec() {
  double x1 = MIN_PRESSURE;
  double f1 = 0.0;
  double x2 = MIN_PRESSURE + PDD_DELTA;
  double f2 =
      std::pow(((x2 - MIN_PRESSURE) / (NORMAL_PRESSURE - MIN_PRESSURE)), 0.5);
  double df1 = PDD_SLOPE;
  double df2 =
      0.5 *
      std::pow(((x2 - MIN_PRESSURE) / (NORMAL_PRESSURE - MIN_PRESSURE)),
               (-0.5)) *
      1.0 / (NORMAL_PRESSURE - MIN_PRESSURE);
  return compute_poly_coefficients({x1, x2}, {f1, f2}, {df1, df2});
}

Eigen::Vector4d pipenetwork::Curves::compute_pdd2_poly_vec() {
  double x1 = NORMAL_PRESSURE - PDD_DELTA;
  double f1 =
      std::pow(((x1 - MIN_PRESSURE) / (NORMAL_PRESSURE - MIN_PRESSURE)), 0.5);
  double x2 = NORMAL_PRESSURE;
  double f2 = 1;
  double df1 =
      0.5 *
      std::pow(((x1 - MIN_PRESSURE) / (NORMAL_PRESSURE - MIN_PRESSURE)),
               (-0.5)) *
      1.0 / (NORMAL_PRESSURE - MIN_PRESSURE);
  double df2 = PDD_SLOPE;
  return compute_poly_coefficients({x1, x2}, {f1, f2}, {df1, df2});
}

Eigen::Vector4d pipenetwork::Curves::compute_leak_poly_vec(double leak_area) {
  double x1 = 0.0;
  double f1 = 0.0;
  double df1 = 1.0e-11;
  double x2 = 1e-4;
  double f2 = LEAK_COEFF * leak_area * std::pow((2 * G * x2), 0.5);
  double df2 = 0.5 * LEAK_COEFF * leak_area * std::pow((2 * G), 0.5) *
               std::pow(x2, -0.5);

  return compute_poly_coefficients({x1, x2}, {f1, f2}, {df1, df2});
}
