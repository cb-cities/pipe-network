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
        if ((current_key != "init") && (info.size() > 0)) {
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


void pipenetwork::Input::construct_node_coord() {
  std::string buf;
  Index id, mesh_id = 0;
  double x_coord, y_coord;

  for (auto const& line : sections_.at("[COORDINATES]")) {
    // skip keys entries
    if (line[0] == '[' || line[0] == ';') continue;

    std::istringstream iss(line);
    // parsing lines
    if (iss >> id >> x_coord >> y_coord) {
      Eigen::Vector3d coord{x_coord, y_coord, 0};
      node_id_map_[id] = mesh_id;
      node_coords_.emplace_back(std::make_pair(mesh_id, coord));
      ++mesh_id;
    }
  }
}

void pipenetwork::Input::construct_node_elevation() {

  // get elevation for junctions
  auto junction_elevation = parse_node_line("[JUNCTIONS]", "elevation");
  njunction_ = junction_elevation.size();

  node_elevations_.insert(std::end(node_elevations_),
                          std::begin(junction_elevation),
                          std::end(junction_elevation));

  // get elevation for reservoirs
  auto reservoir_elevation = parse_node_line("[RESERVOIRS]", "elevation");
  nresvoir_ = reservoir_elevation.size();
  node_elevations_.insert(std::end(node_elevations_),
                          std::begin(reservoir_elevation),
                          std::end(reservoir_elevation));

  // get elevation for reservoirs
  auto tank_elevation = parse_node_line("[TANKS]", "elevation");
  ntank_ = reservoir_elevation.size();
  node_elevations_.insert(std::end(node_elevations_),
                          std::begin(tank_elevation), std::end(tank_elevation));
}

void pipenetwork::Input::construct_node_demand() {

  // get demand for junctions
  auto junction_demand =
      parse_node_line("[JUNCTIONS]", "demand");  // parse junction
  node_demands_.insert(std::end(node_demands_), std::begin(junction_demand),
                       std::end(junction_demand));
  // get demand for reservoir
  auto reservoir_demand =
      parse_node_line("[RESERVOIRS]", "demand");  // parse junction
  node_demands_.insert(std::end(node_demands_), std::begin(reservoir_demand),
                       std::end(reservoir_demand));
  // get demand for tank
  auto tank_demand = parse_node_line("[TANKS]", "demand");  // parse junction
  node_demands_.insert(std::end(node_demands_), std::begin(tank_demand),
                       std::end(tank_demand));
}

std::vector<std::pair<Index, double>> pipenetwork::Input::parse_node_line(
    std::string section_name, std::string mode) const {
  // check if the section name is right.
  if (!sections_.count(section_name)) {
    throw std::invalid_argument(
        "Parsing elevation faild! Invalid section name");
  }

  double elevation, demand;
  Index id, mesh_id;

  std::vector<std::pair<Index, double>> node_info;
  // get elevation for junctions
  for (auto const& line : sections_.at(section_name)) {
    // skip keys entries
    if (line[0] == '[' || line[0] == ';') continue;

    std::istringstream iss(line);
    // parsing lines for demand
    if (mode == "demand") {
      if (section_name == "[JUNCTIONS]") {
        if (iss >> id >> elevation >> demand) {
          mesh_id = node_id_map_.at(id);
          node_info.emplace_back(std::make_pair(mesh_id, to_si (demand,"demand")));
        }

      } else {
        if (iss >> id) {
          mesh_id = node_id_map_.at(id);
          node_info.emplace_back(std::make_pair(mesh_id, -1e-3));  // -1 for
          // reservoir/tank
        }
      }
    }
    // parsing lines for elevation
    else {
      if (iss >> id >> elevation) {
        mesh_id = node_id_map_.at(id);
        node_info.emplace_back(std::make_pair(mesh_id, to_si (elevation,"elevation")));

      }
    }
  }
  return node_info;
}

void pipenetwork::Input::construct_pipe_info() {
  Index pid, nid1, nid2, mesh_id = 0;
  double length, diameter, roughness, loss;
  std::string status;

  // get elevation for junctions
  for (auto const& line : sections_.at("[PIPES]")) {
    // skip keys entries
    if (line[0] == '[' || line[0] == ';') continue;

    std::istringstream iss(line);
    // parsing lines for pipe
    if (iss >> pid >> nid1 >> nid2 >> length >> diameter >> roughness >> loss >>
        status) {

      pipe_id_map_[pid] = mesh_id;
      nodeids_.emplace_back(
          std::make_pair(node_id_map_.at(nid1), node_id_map_.at(nid2)));
      diameter_.emplace_back(to_si (diameter,"diameter"));
      length_.emplace_back(to_si (length,"length"));
      roughness_.emplace_back(roughness);
      if (status == "Open") {
        pipe_status_.emplace_back(true);
      } else {
        pipe_status_.emplace_back(false);
      }

      ++mesh_id;
    }
  }
}


double pipenetwork::to_si (double val, std::string mode) {
    if (mode=="elevation"||mode=="length"){
        return val*0.3048;//ft to meter
    }
    else if (mode == "demand"){
        return val*6.30901964e-05;//GPM to si
    }
    else if (mode == "diameter"){
        return val*0.0254;//inch to meter
    }

}
