#include <node.h>

Eigen::VectorXd pipenetwork::Node::compute_poly_coefficients(
    double x1, double x2, double f1, double f2, double df1, double df2) {
  Eigen::VectorXd ret(4);
  double a = (2 * (f1 - f2) - (x1 - x2) * (df2 + df1)) /
             (std::pow(x2, 3) - std::pow(x1, 3) + 3 * x1 * x2 * (x1 - x2));
  double b = (df1 - df2 + 3 * ((std::pow(x2, 2) - std::pow(x1, 2)) * a)) /
             (2 * (x1 - x2));
  double c = df2 - 3 * std::pow(x2, 2) * a - 2 * x2 * b;
  double d = f2 - std::pow(x2, 3) * a - std::pow(x2, 2) * b - x2 * c;
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
    double df2 = 0.5*std::pow ((x2-minimum_pressure_)/(normal_pressure_-minimum_pressure_),-0.5)
          *(1/(normal_pressure_-minimum_pressure_));

  return compute_poly_coefficients (x1,x2,f1,f2,df1,df2);
}

Eigen::VectorXd pipenetwork::Node::get_pdd_poly_coef_2() {
    double x1 = normal_pressure_-pdd_smoothing_delta_;
    double f1 = std::pow((x1-minimum_pressure_)/(normal_pressure_-minimum_pressure_),0.5);
    double x2 = normal_pressure_;
    double f2 = 1;
    double df1 = 0.5*std::pow((x1-minimum_pressure_)/(normal_pressure_-minimum_pressure_),-0.5)
            *(1/(normal_pressure_-minimum_pressure_));
    double df2 = pdd_slope_;
  return compute_poly_coefficients (x1,x2,f1,f2,df1,df2);
}
