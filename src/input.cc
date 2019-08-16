#include "input.h"

void pipenetwork::Input::parse_sections() {
  std::string line;
  std::string buf, current_key{"init"};
  std::vector<std::string> info;

  std::ifstream infile(filename_);
  //! parsing line one by one
  bool parsing_flag = false;
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    while (iss >> buf) {
      if (buf[0] == '[') {
        // check the condition
        if ((current_key != "init") && (!info.empty())) {
          // put parsed lines into map and clear the buffer vector
          sections_[current_key] = info;
          info.clear();
          parsing_flag = false;
        }
        // update the keys
        bool is_in = section_keys_.find(buf) != section_keys_.end();
        if (is_in) {
          current_key = buf;
          parsing_flag = true;
        }
      }
    }
    // parse the line
    if (parsing_flag) {
      info.emplace_back(line);
    }
  }
}

// void pipenetwork::Input::construct_node_coord() {
//  std::string buf;
//  Index id, mesh_id = 0;
//  double x_coord, y_coord;
//
//  for (auto const& line : sections_.at("[COORDINATES]")) {
//    // skip keys entries
//    if (line[0] == '[' || line[0] == ';') continue;
//
//    std::istringstream iss(line);
//    // parsing lines
//    if (iss >> id >> x_coord >> y_coord) {
//      Eigen::Vector3d coord{x_coord, y_coord, 0};
//      node_coords_.emplace_back(std::make_pair(mesh_id, coord));
//      ++mesh_id;
//    }
//  }
//}

void pipenetwork::Input::construct_node_elevation_head() {

  // get elevation for junctions
  auto junction_info = parse_node_line("[JUNCTIONS]", "elevation");
  junction_ids_ = junction_info.first;
  junction_elevations_ = junction_info.second;

  // get head for reservoirs
  auto reservoir_info = parse_node_line("[RESERVOIRS]", "head");
  reservoir_ids_ = reservoir_info.first;
  reservoir_heads_ = reservoir_info.second;
}

void pipenetwork::Input::construct_node_demand() {
  // get demand for junctions
  auto junction_demand_info =
      parse_node_line("[JUNCTIONS]", "demand");  // parse junction
  junction_ids_ = junction_demand_info.first;
  junction_demands_ = junction_demand_info.second;
}

std::pair<std::vector<std::string>, std::vector<double>>
    pipenetwork::Input::parse_node_line(const std::string& section_name,
                                        const std::string& mode) const {
  // check if the section name is right.

  if (!sections_.count(section_name)) {
    throw std::invalid_argument("Parsing faild! Invalid section name");
  }

  double elevation_head, demand;
  std::string id;

  std::vector<std::string> node_ids;
  std::vector<double> node_info;
  // get elevation for junctions
  for (auto const& line : sections_.at(section_name)) {
    // skip keys entries
    if (line[0] == '[' || line[0] == ';') continue;

    std::istringstream iss(line);
    // parsing lines for demand
    if (mode == "demand") {
      if (iss >> id >> elevation_head >> demand) {
        node_info.emplace_back(to_si(demand, "demand"));
        node_ids.emplace_back(id);
      }
    }
    // parsing lines for elevation/head
    else {
      if (iss >> id >> elevation_head) {
        node_info.emplace_back(to_si(elevation_head, "elevation"));
        node_ids.emplace_back(id);
      }
    }
  }

  return std::make_pair(node_ids, node_info);
}

void pipenetwork::Input::construct_node_info() {
  construct_node_elevation_head();
  construct_node_demand();

  for (int i = 0; i < junction_elevations_.size(); ++i) {
    pipenetwork::Junction_prop junc_prop;
    junc_prop.id = junction_ids_[i];
    junc_prop.elevation = junction_elevations_[i];
    junc_prop.demand = junction_demands_[i];

    junc_props_.emplace_back(junc_prop);
  }

  for (int i = 0; i < reservoir_ids_.size(); ++i) {
    pipenetwork::Reservoir_prop res_prop;
    res_prop.id = reservoir_ids_[i];
    res_prop.head = reservoir_heads_[i];

    res_props_.emplace_back(res_prop);
  }
}

void pipenetwork::Input::construct_pipe_info() {
  Index mesh_id = 0;
  double length, diameter, roughness, loss;
  std::string pid, nid1, nid2, status;

  // get pipe information
  for (auto const& line : sections_.at("[PIPES]")) {
    // skip keys entries
    if (line[0] == '[' || line[0] == ';') continue;

    std::istringstream iss(line);
    // parsing lines for pipe
    if (iss >> pid >> nid1 >> nid2 >> length >> diameter >> roughness >> loss >>
        status) {
      std::transform(status.begin(), status.end(), status.begin(), ::toupper);
      pipe_ids_.emplace_back(pid);
      nodeids_.emplace_back(std::make_pair(nid1, nid2));
      diameter_.emplace_back(to_si(diameter, "diameter"));
      length_.emplace_back(to_si(length, "length"));
      roughness_.emplace_back(roughness);

      if (status == "OPEN") {
        pipe_status_.emplace_back(OPEN);
      } else {
        pipe_status_.emplace_back(CLOSED);
      }
      ++mesh_id;
    }
  }

  for (int i = 0; i < pipe_ids_.size(); ++i) {
    pipenetwork::Pipe_prop pipe_prop;
    pipe_prop.id = pipe_ids_[i];
    pipe_prop.length = length_[i];
    pipe_prop.diameter = diameter_[i];
    pipe_prop.roughness = roughness_[i];
    pipe_prop.node1_id = nodeids_[i].first;
    pipe_prop.node2_id = nodeids_[i].second;
    pipe_prop.status = pipe_status_[i];
    pipe_props_.emplace_back(pipe_prop);
  }
}

void pipenetwork::Input::construct_pump_info() {
  // get pump information
  std::string pid, nid1, nid2, type, curve_name;
  double speed{1.0}, power{50};
  for (auto const& line : sections_.at("[PUMPS]")) {
    Pump_prop pump_prop;
    // skip keys entries
    if (line[0] == '[' || line[0] == ';') continue;

    std::istringstream iss(line);

    if (iss >> pid >> nid1 >> nid2 >> type) {
      std::transform(type.begin(), type.end(), type.begin(), ::toupper);
      if (type == "HEAD") {
        iss >> curve_name;
        pump_prop.curve_name = curves_info_->pump_str_int(curve_name);
        pump_prop.id = pid;
        pump_prop.node1_id = nid1;
        pump_prop.node2_id = nid2;
        pump_prop.pump_type = HEADPUMP;
      }
      if (type == "POWER") {
        iss >> power;
        pump_prop.id = pid;
        pump_prop.node1_id = nid1;
        pump_prop.node2_id = nid2;
        pump_prop.pump_type = POWERPUMP;
        pump_prop.power = to_si(power, "power");
      }
      pump_props_.emplace_back(pump_prop);
    }
  }
}

void pipenetwork::Input::construct_valve_info() {
  // get valve information
  std::string pid, nid1, nid2, type;
  double diameter, setting, minorloss{0};
  for (auto const& line : sections_.at("[VALVES]")) {
    Valve_prop valve_prop;
    // skip keys entries
    if (line[0] == '[' || line[0] == ';') continue;
    std::istringstream iss(line);
    if (iss >> pid >> nid1 >> nid2 >> diameter >> type >> setting) {
      std::transform(type.begin(), type.end(), type.begin(), ::toupper);
      valve_prop.id = pid;
      valve_prop.node1_id = nid1;
      valve_prop.node2_id = nid2;
      valve_prop.diameter = to_si(diameter, "diameter");

      if (type == "PRV") {
        valve_prop.valve_type = PRVALVE;
        valve_prop.setting = to_si(setting, "pressure");
      }
      if (type == "FCV") {
        valve_prop.valve_type = FCVALVE;
        valve_prop.setting = to_si(setting, "flow");
      }
      if (type == "TCV") {
        valve_prop.valve_type = TCVALVE;
        valve_prop.setting = setting;
      }
      if (iss >> minorloss) {
        valve_prop.minor_loss = minorloss;
      }
      valve_props_.emplace_back(valve_prop);
    }
  }
}

void pipenetwork::Input::construct_curve_info() {
  // get curve information
  std::string curve_id;
  double x, y;
  std::string id_buff;
  std::vector<std::pair<double, double>> curve_point;
  std::vector<Pump_curve_prop> head_pump_props;

  for (auto const& line : sections_.at("[CURVES]")) {
    // skip keys entries
    if (line[0] == '[' || line[0] == ';') continue;
    std::istringstream iss(line);

    if (iss >> curve_id >> x >> y) {
      if (!curve_point.empty() && curve_id != id_buff) {
        Pump_curve_prop pump_curve(id_buff, curve_point);
        head_pump_props.emplace_back(pump_curve);
        curve_point.clear();
      }
      curve_point.emplace_back(
          std::make_pair(to_si(x, "flow"), to_si(y, "head")));
      id_buff = curve_id;
    }
  }
  if (!curve_point.empty()) {
    Pump_curve_prop pump_curve(curve_id, curve_point);
    head_pump_props.emplace_back(pump_curve);
    // add pump info to curve object
    curves_info_->add_pump_curves(head_pump_props);
  }
}

double pipenetwork::to_si(double val, const std::string& mode) {
  if (mode == "elevation" || mode == "length" || mode == "head") {
    return val * 0.3048;  // ft to meter
  } else if (mode == "demand" || mode == "flow") {
    return val * 6.30901964e-05;  // GPM to si
  } else if (mode == "diameter") {
    return val * 0.0254;  // inch to meter
  } else if (mode == "power") {
    return val * 745.699872;  // hp to W (Nm/s)
  } else if (mode == "pressure") {
    return val * (0.3048 / 0.4333);  // psi * (m/ft / psi/ft)
  } else {
    throw std::runtime_error("Mode not recognized!");
  }
}
