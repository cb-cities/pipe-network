#include "io.h"

void pipenetwork::IO::read_inp(const std::string& filename) {
  filename_ = filename;
  if (!IO_utils::file_exists(filename)) {
    throw std::runtime_error(
        "Input file does not exist, please check the path!");
  }
  try {
    parse_sections();
    construct_node_info();
    construct_pipe_info();
    construct_curve_info();
    construct_pump_info();
    construct_valve_info();
  } catch (std::exception& e) {
    std::cerr << "Failed to read input file, error message " << e.what()
              << std::endl;
    std::abort();
  }
}

void pipenetwork::IO::create_synthetic_net(pipenetwork::Index n) {
  auto junction_nodes = construct_synthesis_junctions(n);
  construct_synthesis_pipes(junction_nodes);
  create_sources(n);
}

void pipenetwork::IO::parse_sections() {
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

std::pair<std::vector<std::string>, std::vector<double>>
    pipenetwork::IO::parse_node_line(const std::string& section_name,
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
        node_info.emplace_back(IO_utils::to_si(demand, "demand"));
        node_ids.emplace_back(id);
      }
    }
    // parsing lines for elevation/head
    else {
      if (iss >> id >> elevation_head) {
        node_info.emplace_back(IO_utils::to_si(elevation_head, "elevation"));
        node_ids.emplace_back(id);
      }
    }
  }

  return std::make_pair(node_ids, node_info);
}

void pipenetwork::IO::parse_leak_info() {
  std::string nid;
  double diameter;
  for (auto const& line : sections_.at("[LEAKS]")) {
    // skip keys entries
    if (line[0] == '[' || line[0] == ';') continue;

    std::istringstream iss(line);
    if (iss >> nid >> diameter) {
      nid2ldia_[nid] = diameter;
    }
  }
}

void pipenetwork::IO::construct_node_info() {
  // get elevation for junctions
  auto junction_info = parse_node_line("[JUNCTIONS]", "elevation");
  auto junction_ids = junction_info.first;
  auto junction_elevations = junction_info.second;
  // get demand for junctions
  auto junction_demand_info =
      parse_node_line("[JUNCTIONS]", "demand");  // parse junction
  auto junction_demands = junction_demand_info.second;

  // get head for reservoirs
  auto reservoir_info = parse_node_line("[RESERVOIRS]", "head");
  auto reservoir_ids = reservoir_info.first;
  auto reservoir_heads = reservoir_info.second;

  // get leak information
  bool leak = sections_.find("[LEAKS]") != sections_.end();
  if (leak) {
    parse_leak_info();
  }

  // construct property vectors
  for (int i = 0; i < junction_elevations.size(); ++i) {
    pipenetwork::JunctionProp junc_prop;
    junc_prop.name = junction_ids[i];
    junc_prop.elevation = junction_elevations[i];
    junc_prop.demand = junction_demands[i];
    bool leak_node = nid2ldia_.find(junction_ids[i]) != nid2ldia_.end();
    if (leak_node) {
      junc_prop.leak_diameter =
          IO_utils::to_si(nid2ldia_.at(junction_ids[i]), "diameter");
    }
    junc_props_.emplace_back(junc_prop);
  }

  for (int i = 0; i < reservoir_ids.size(); ++i) {
    pipenetwork::ReservoirProp res_prop;
    res_prop.name = reservoir_ids[i];
    res_prop.head = reservoir_heads[i];
    res_props_.emplace_back(res_prop);
  }
}

void pipenetwork::IO::construct_pipe_info() {
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

      pipenetwork::PipeProp pipe_prop;
      std::transform(status.begin(), status.end(), status.begin(), ::toupper);
      pipe_prop.name = pid;
      pipe_prop.length = IO_utils::to_si(length, "length");
      pipe_prop.diameter = IO_utils::to_si(diameter, "diameter");
      pipe_prop.roughness = roughness;
      pipe_prop.node1_name = nid1;
      pipe_prop.node2_name = nid2;

      if (status == "OPEN") {
        pipe_prop.status = LinkStatus::OPEN;
      } else {
        pipe_prop.status = LinkStatus::CLOSED;
      }
      pipe_props_.emplace_back(pipe_prop);
    }
  }
}

void pipenetwork::IO::construct_pump_info() {
  // get pump information
  std::string pid, nid1, nid2, type, curve_name;
  double power{50};
  for (auto const& line : sections_.at("[PUMPS]")) {
    PumpProp pump_prop;
    // skip keys entries
    if (line[0] == '[' || line[0] == ';') continue;

    std::istringstream iss(line);

    if (iss >> pid >> nid1 >> nid2 >> type) {
      std::transform(type.begin(), type.end(), type.begin(), ::toupper);
      if (type == "HEAD") {
        iss >> curve_name;
        pump_prop.curve_id = curves_info_->pump_str_int(curve_name);
        pump_prop.name = pid;
        pump_prop.node1_name = nid1;
        pump_prop.node2_name = nid2;
        pump_prop.type = PumpType::HEADPUMP;
      }
      if (type == "POWER") {
        iss >> power;
        pump_prop.name = pid;
        pump_prop.node1_name = nid1;
        pump_prop.node2_name = nid2;
        pump_prop.type = PumpType::POWERPUMP;
        pump_prop.power = IO_utils::to_si(power, "power");
      }
      pump_props_.emplace_back(pump_prop);
    }
  }
}

void pipenetwork::IO::construct_valve_info() {
  // get valve information
  std::string pid, nid1, nid2, type;
  double diameter, setting, minorloss{0};
  for (auto const& line : sections_.at("[VALVES]")) {
    ValveProp valve_prop;
    // skip keys entries
    if (line[0] == '[' || line[0] == ';') continue;
    std::istringstream iss(line);
    if (iss >> pid >> nid1 >> nid2 >> diameter >> type >> setting) {
      std::transform(type.begin(), type.end(), type.begin(), ::toupper);
      valve_prop.name = pid;
      valve_prop.node1_name = nid1;
      valve_prop.node2_name = nid2;
      valve_prop.diameter = IO_utils::to_si(diameter, "diameter");
      valve_prop.status = LinkStatus::ACTIVE;

      if (type == "PRV") {
        valve_prop.type = ValveType::PRVALVE;
        valve_prop.setting = IO_utils::to_si(setting, "pressure");
      } else if (type == "FCV") {
        valve_prop.type = ValveType::FCVALVE;
        valve_prop.setting = IO_utils::to_si(setting, "flow");
      } else if (type == "TCV") {
        valve_prop.type = ValveType::TCVALVE;
        valve_prop.setting = setting;
      }
      if (iss >> minorloss) {
        valve_prop.minor_loss_coeff = minorloss;
      }
      valve_props_.emplace_back(valve_prop);
    }
  }
}

void pipenetwork::IO::construct_curve_info() {
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
      curve_point.emplace_back(std::make_pair(IO_utils::to_si(x, "flow"),
                                              IO_utils::to_si(y, "head")));
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

std::vector<std::vector<std::string>>
    pipenetwork::IO::construct_synthesis_junctions(int n) {
  // buffer vector
  std::vector<std::vector<std::string>> junction_names(n);
  std::vector<std::string> buffer;
  std::string junction_name;
  pipenetwork::JunctionProp junc_prop;
  double elevation, demand;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      // create junction name
      junction_name = "J-" + std::to_string(i * n + j);
      // elevation and demand
      elevation = IO_utils::rand_number(400, 900);
      demand = IO_utils::rand_number(0, 5);
      // add to junction array
      junc_prop.name = junction_name;
      junc_prop.elevation = IO_utils::to_si(elevation, "elevation");
      junc_prop.demand = IO_utils::to_si(demand, "demand");
      junc_props_.emplace_back(junc_prop);
      buffer.emplace_back(junction_name);
    }
    junction_names[i] = buffer;
    buffer.clear();
  }
  return junction_names;
}

void pipenetwork::IO::construct_synthesis_pipes(
    const std::vector<std::vector<std::string>>& junction_names) {
  int n = junction_names.size();
  for (int i = 0; i < n; ++i) {
    std::vector<std::string> junction_name_col = junction_names.at(i);
    create_vertical_pipes(junction_name_col, i);
  }
  for (int i = 1; i < n; ++i) {
    std::vector<std::string> left_junc = junction_names.at(i - 1);
    std::vector<std::string> right_junc = junction_names.at(i);
    create_horizontal_pipes(left_junc, right_junc, i);
  }
}

void pipenetwork::IO::create_vertical_pipes(
    const std::vector<std::string>& junction_names, int col_num, bool rand) {
  int n = junction_names.size();
  pipenetwork::PipeProp pipe_prop;
  std::string upper_junc_id, lower_junc_id, pipe_name;
  for (int i = 1; i < n; ++i) {
    if (rand) {
      if (IO_utils::rand_number(0, 1) < 0.5) {
        continue;
      }
    }

    lower_junc_id = junction_names.at(i - 1);
    upper_junc_id = junction_names.at(i);
    pipe_name = "P-V-" + std::to_string(col_num * n + i);
    pipe_prop.name = pipe_name;
    pipe_prop.length =
        IO_utils::to_si(IO_utils::rand_number(200, 800), "length");
    pipe_prop.diameter =
        IO_utils::to_si(IO_utils::rand_number(3, 15), "diameter");
    pipe_prop.roughness = 155;
    pipe_prop.node1_name = lower_junc_id;
    pipe_prop.node2_name = upper_junc_id;
    pipe_props_.emplace_back(pipe_prop);
  }
}

void pipenetwork::IO::create_horizontal_pipes(
    const std::vector<std::string>& l_junc,
    const std::vector<std::string>& r_junc, int col_num, bool rand) {
  int n = l_junc.size();
  pipenetwork::PipeProp pipe_prop;
  std::string left_junc_id, right_junc_id, pipe_name;
  for (int i = 0; i < n; ++i) {
    if (rand) {
      if (IO_utils::rand_number(0, 1) < 0.5) {
        continue;
      }
    }
    left_junc_id = l_junc.at(i);
    right_junc_id = r_junc.at(i);
    pipe_name = "P-H-" + std::to_string(col_num * n + i);
    pipe_prop.name = pipe_name;
    pipe_prop.length =
        IO_utils::to_si(IO_utils::rand_number(200, 800), "length");
    pipe_prop.diameter =
        IO_utils::to_si(IO_utils::rand_number(3, 15), "diameter");
    pipe_prop.roughness = 155;
    pipe_prop.node1_name = left_junc_id;
    pipe_prop.node2_name = right_junc_id;
    pipe_props_.emplace_back(pipe_prop);
  }
}

void pipenetwork::IO::create_sources(int n) {
  pipenetwork::ReservoirProp src_prop;
  pipenetwork::PipeProp pipe_prop;
  std::string src_name, pipe_name;

  for (int i = 0; i < 5; ++i) {
    src_name = "SRC-" + std::to_string(i);
    src_prop.name = src_name;
    src_prop.head =
        IO_utils::to_si(IO_utils::rand_number(500, 1000), "elevation");

    pipe_name = "P-S-" + std::to_string(i);
    pipe_prop.name = pipe_name;
    pipe_prop.length = IO_utils::to_si(30, "length");
    pipe_prop.diameter = IO_utils::to_si(30, "diameter");
    pipe_prop.roughness = 155;
    pipe_prop.node1_name = src_name;
    pipe_prop.node2_name = "J-" + std::to_string(i * n);

    res_props_.emplace_back(src_prop);
    pipe_props_.emplace_back(pipe_prop);
  }
}

void pipenetwork::IO::save_mesh_inp(
    const std::shared_ptr<pipenetwork::Mesh>& mesh,
    const std::string& output_path) {
  auto out = IO_utils::Output(mesh);
  out.save_mesh_inp(output_path);
}

void pipenetwork::IO::save_sim_result(
    const std::shared_ptr<pipenetwork::Mesh>& mesh,
    const std::string& output_path) {
  auto out = IO_utils::Output(mesh);
  out.save_sim_result(output_path);
}
