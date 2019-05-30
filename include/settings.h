#ifndef PIPE_NETWORK_SETTING_H_
#define PIPE_NETWORK_SETTING_H_

//! define alias of type for index
using Index = unsigned long long;

enum Node_type { junction, tank, reservoir };
enum Pipe_type { pipe, valve, pump };

namespace pipenetwork {
//! define gravitational acceleration in x, y, z direction
const Eigen::Vector3d Gravity(0.00, 0.00, 9.81);
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_SETTING_H_
