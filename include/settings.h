#ifndef PIPE_NETWORK_SETTING_H_
#define PIPE_NETWORK_SETTING_H_

//! define alias of type for index
using Index = unsigned long long;


//! node settings
enum Node_type { junction_type, tank_type, reservoir_type };
const double PI{3.1415926};
const double G{9.81};


//! link settings
enum Link_type { pipe_type, valve_type, pump_type };
enum Pipe_status {open,closed};

#endif  // PIPE_NETWORK_SETTING_H_
