#ifndef PIPE_NETWORK_SETTING_H_
#define PIPE_NETWORK_SETTING_H_
#include <cmath>
//! define alias of type for index
using Index = unsigned long long;

//! node settings
enum Node_type { JUNCTION, TANK, RESERVOIR };
const double PI{M_PI};
const double G{9.81};
const double HW_COEFF{10.666829500036352};
const double LEAK_COEFF{0.75};
//! pdd settings
const double MIN_PRESSURE{-1e-12};
const double NORMAL_PRESSURE{50};
const Eigen::Vector4d PDD_POLY_VEC1{-18.749999999749996, 6.2499999999,
                                    1.000000082740371e-11,
                                    -4.440892098516782e-17};
const Eigen::Vector4d PDD_POLY_VEC2{-0.6249920885505783, 37.249212823040864,
                                    -739.9780066609305, 4900.811712406892};
const double PDD_DELTA{0.2};
const double PDD_SLOPE{1e-12};

//! link settings
enum Link_type { PIPE, VALVE, PUMP };
enum Pipe_status { OPEN, CLOSED };

//! constant for modified hazen-williams formula
const double HW_Q1{0.0002};
const double HW_Q2{0.0004};
const double HW_M{0.001};
const Eigen::Vector4d HW_POLY_VEC{6619.952473405493, -2.562247355522429,
                                  0.0012305046454003125,
                                  3.4293453535907055e-09};

#endif  // PIPE_NETWORK_SETTING_H_
