#include "catch.hpp"
#include <iostream>
#include <string>

#include "curves.h"
// Check node class
TEST_CASE("Check Curves", "[Curve]") {
  // Tolerance
  const double tolerance = 1.e-8;
  std::string head_pump_name1 = "1";
  std::string head_pump_name2 = "2";

  std::vector<std::pair<double, double>> pump_curve_points1 = {
      {0.0, 31.6992},
      {0.1261803928, 28.041600000000003},
      {0.2523607856, 19.2024}};
  std::vector<std::pair<double, double>> pump_curve_points2 = {
      {0.1261803928, 30.48}};
  pipenetwork::Pump_curve_prop p_curve1(head_pump_name1, pump_curve_points1);
  pipenetwork::Pump_curve_prop p_curve2(head_pump_name2, pump_curve_points2);
  pipenetwork::Curves curve_list;

  SECTION("Check head pump curve property for three points input") {

    auto head_curve_coefficients = p_curve1.head_curve_coefficients;

    // check coefficients
    REQUIRE(head_curve_coefficients[0] == Approx(31.6992).epsilon(tolerance));
    REQUIRE(head_curve_coefficients[1] ==
            Approx(143.47246994481017).epsilon(tolerance));
    REQUIRE(head_curve_coefficients[2] ==
            Approx(1.7725895038969284).epsilon(tolerance));
  }

  SECTION("Check head pump curve property for one point input") {

    auto head_curve_coefficients = p_curve2.head_curve_coefficients;

    // check coefficients
    REQUIRE(head_curve_coefficients[0] == Approx(40.64).epsilon(tolerance));
    REQUIRE(head_curve_coefficients[1] ==
            Approx(638.1311689716313).epsilon(tolerance));
    REQUIRE(head_curve_coefficients[2] == Approx(2).epsilon(tolerance));
  }

  SECTION("Check linear approximation coefficients") {
    auto line = p_curve1.line_param;
    REQUIRE(line[0] == Approx(4.4542219126875336e-18).epsilon(tolerance));
    REQUIRE(line[1] == Approx(31.6992).epsilon(tolerance));
  }
  SECTION("Check add pump properties") {
    std::vector<pipenetwork::Pump_curve_prop> p_curve_vec = {p_curve1,
                                                             p_curve2};
    curve_list.add_pump_curves(p_curve_vec);
    REQUIRE(curve_list.pump_curves()["1"].line_param[1] ==
            Approx(31.6992).epsilon(tolerance));
    REQUIRE(curve_list.poly_coeffs()["HW_POLY_VEC"][0] ==
            Approx(6619.952473405493).epsilon(tolerance));
    REQUIRE(curve_list.poly_coeffs()["HW_POLY_VEC"][1] ==
            Approx(-2.562247355522429).epsilon(tolerance));
    REQUIRE(curve_list.poly_coeffs()["HW_POLY_VEC"][2] ==
            Approx(0.0012305046454003125).epsilon(tolerance));
    REQUIRE(curve_list.poly_coeffs()["HW_POLY_VEC"][3] ==
            Approx(3.4293453535907055e-09).epsilon(tolerance));

    //    REQUIRE(curve_list.poly_coeffs()["PDD_POLY_VEC1"][0] ==
    //            Approx(-18.749999999749996).epsilon(tolerance));
    //    REQUIRE(curve_list.poly_coeffs()["PDD_POLY_VEC1"][1] ==
    //            Approx(6.2456331391).epsilon(tolerance));
    //    REQUIRE(curve_list.poly_coeffs()["PDD_POLY_VEC1"][2] ==
    //            Approx(1.000000082740371e-11).epsilon(tolerance));
    //    REQUIRE(curve_list.poly_coeffs()["PDD_POLY_VEC1"][3] ==
    //            Approx(-4.440892098516782e-17).epsilon(tolerance));
    //
    //    REQUIRE(curve_list.poly_coeffs()["PDD_POLY_VEC2"][0] ==
    //            Approx(-0.6249920885505783).epsilon(tolerance));
    //    REQUIRE(curve_list.poly_coeffs()["PDD_POLY_VEC2"][1] ==
    //            Approx(37.249212823040864).epsilon(tolerance));
    //    REQUIRE(curve_list.poly_coeffs()["PDD_POLY_VEC2"][2] ==
    //            Approx(-739.9780066609305).epsilon(tolerance));
    //    REQUIRE(curve_list.poly_coeffs()["PDD_POLY_VEC2"][3] ==
    //            Approx(4900.811712406892).epsilon(tolerance));
  }
}