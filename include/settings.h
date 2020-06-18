#ifndef PIPE_NETWORK_SETTING_H_
#define PIPE_NETWORK_SETTING_H_
#include <cmath>
namespace pipenetwork {
//! define alias of type for index
using Index = unsigned long;

//! node settings
const double PI{M_PI};
const double G{9.81};
const double HW_COEFF{10.666829500036352};
const double LEAK_COEFF{0.75};
const double INIT_FLOWRATE{1E-3};

//! pressure demand driven simulation (PDD) settings
const double MIN_PRESSURE{1e-4};
const double NORMAL_PRESSURE{20};
const double PDD_DELTA{0.2};
const double PDD_SLOPE{1e-11};

//! link settings
enum class PumpType { POWERPUMP, HEADPUMP };
enum class ValveType { PRVALVE, FCVALVE, TCVALVE };
enum class LinkStatus { OPEN, CLOSED, ACTIVE };
const double PUMP_POWER{50};
const double PUMP_SPEED{1};
const double MINOR_LOSS_COEFF{0};
const double VALVE_DIAMETER{0.3048};

//! constant for modified hazen-williams formula
const double HW_Q1{0.0002};
const double HW_Q2{0.0004};
const double HW_M{0.001};

//! pump settings
const double PUMP_M = -0.00000000001;
const double PUMP_Q1 = 0.0;
const double PUMP_Q2 = 1.0e-8;
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_SETTING_H_
