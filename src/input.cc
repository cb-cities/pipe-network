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

std::pair<std::vector<Index>, std::vector<double>>
    pipenetwork::Input::parse_node_line(const std::string& section_name,
                                        const std::string& mode) const {
  // check if the section name is right.
  if (!sections_.count(section_name)) {
    throw std::invalid_argument(
        "Parsing elevation faild! Invalid section name");
  }

  double elevation_head, demand;
  Index id;

  std::vector<Index> node_ids;
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

void pipenetwork::Input::construct_pipe_info() {
  Index pid, nid1, nid2, mesh_id = 0;
  double length, diameter, roughness, loss;
  std::string status;

  // get pipe information
  for (auto const& line : sections_.at("[PIPES]")) {
    // skip keys entries
    if (line[0] == '[' || line[0] == ';') continue;

    std::istringstream iss(line);
    // parsing lines for pipe
    if (iss >> pid >> nid1 >> nid2 >> length >> diameter >> roughness >> loss >>
        status) {
      pipe_ids_.emplace_back(pid);
      nodeids_.emplace_back(std::make_pair(nid1, nid2));
      diameter_.emplace_back(to_si(diameter, "diameter"));
      length_.emplace_back(to_si(length, "length"));
      roughness_.emplace_back(roughness);
      if (status == "Open") {
        pipe_status_.emplace_back(OPEN);
      } else {
        pipe_status_.emplace_back(CLOSED);
      }
      ++mesh_id;
    }
  }
}

double pipenetwork::to_si(double val, const std::string& mode) {
  if (mode == "elevation" || mode == "length") {
    return val * 0.3048;  // ft to meter
  } else if (mode == "demand") {
    return val * 6.30901964e-05;  // GPM to si
  } else if (mode == "diameter") {
    return val * 0.0254;  // inch to meter
  } else {
    throw std::runtime_error("Mode not recognized!");
  }
}
