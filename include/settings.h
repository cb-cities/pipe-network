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


//! link settings
enum Link_type { PIPE, VALVE, PUMP };
enum Pipe_status {OPEN,CLOSED};

//! constant for modified hazen-williams formula
const double HW_Q1{0.0002};
const double HW_Q2{0.0004};
const double HW_M{0.001};
const Eigen::Vector4d HW_POLY_VEC{6619.952473405493,-2.562247355522429,0.0012305046454003125,3.4293453535907055e-09};

#endif  // PIPE_NETWORK_SETTING_H_
