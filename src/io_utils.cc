#include "io_utils.h"

double pipenetwork::IO_utils::rand_number(double l, double h) {
  double r = l + static_cast<float>(std::rand()) /
                     (static_cast<float>(RAND_MAX / (h - l)));
  return r;
}

double pipenetwork::IO_utils::to_si(double val, const std::string& mode) {
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

double pipenetwork::IO_utils::from_si(double val, const std::string& mode) {
  if (mode == "elevation" || mode == "length" || mode == "head") {
    return val / 0.3048;  // meter to ft
  } else if (mode == "demand" || mode == "flow") {
    return val / (6.30901964e-05);  // si to GPM
  } else if (mode == "diameter") {
    return val / 0.0254;  // meter to inch
  } else if (mode == "power") {
    return val / 745.699872;  // W (Nm/s) to hp
  } else if (mode == "pressure") {
    return val / (0.3048 / 0.4333);  // psi * (m/ft / psi/ft)
  } else {
    throw std::runtime_error("Mode not recognized!");
  }
}

bool pipenetwork::IO_utils::path_exist(const std::string& pathname) {
  struct stat info;
  if (stat(pathname.c_str(), &info) != 0) {
    return false;
  }
  return info.st_mode & S_IFDIR;
}

void pipenetwork::IO_utils::create_new_folder(const std::string& pathname) {
  if (!path_exist(pathname)) {
    std::cout << "output folder does not exist, creating a new one"
              << std::endl;
    try {
      mkdir(pathname.c_str(), 0777);
    } catch (...) {
      std::cerr << "Failed to create folder, abort";
      abort();
    }
  }
}

void pipenetwork::IO_utils::Output::save_mesh_inp(
    const std::string& output_path) {
  create_new_folder(output_path);
  filepath_ = output_path + mesh_->name() + ".inp";
  std::ofstream outfile{filepath_};
  write_junctions_inp(outfile);
  write_reservoirs_inp(outfile);
  write_pipes_inp(outfile);
}

void pipenetwork::IO_utils::Output::write_junctions_inp(
    std::ofstream& outfile) {
  // write row info
  outfile << "[JUNCTIONS]"
          << "\n";
  outfile << ";ID    "
          << "Elev   "
          << "Demand     "
          << "Pattern     "
          << "\n";
  auto junctions_map = mesh_->nodes()->junctions();
  for (const auto& node : junctions_map) {
    auto junction = node.second;
    auto junction_prop = junction->property();
    outfile << std::setprecision(12) << std::to_string(junction->id()) << "   "
            << from_si(junction_prop.elevation, "elevation") << "   "
            << from_si(junction_prop.demand, "demand") << "   "
            << "        ;"
            << "\n";
  }
}

void pipenetwork::IO_utils::Output::write_reservoirs_inp(
    std::ofstream& outfile) {
  outfile << "\n";
  outfile << "[RESERVOIRS]"
          << "\n";
  outfile << ";ID    "
          << "Head   "
          << "Pattern     "
          << "\n";
  auto res_map = mesh_->nodes()->reservoirs();
  for (const auto& node : res_map) {
    // get node information, determine possible leak nodes and assemble demands
    // heads vector
    auto res = node.second;
    auto res_property = res->property();

    outfile << std::setprecision(12) << std::to_string(res->id()) << "   "
            << from_si(res_property.head, "head") << "   "
            << "        ;"
            << "\n";
  }
}

void pipenetwork::IO_utils::Output::write_pipes_inp(std::ofstream& outfile) {
  outfile << "\n";
  outfile << "[PIPES]"
          << "\n";
  outfile << ";ID    "
          << "Node1   "
          << "Node2   "
          << "Length   "
          << "Diameter   "
          << "Roughness   "
          << "MinorLoss   "
          << "Status   "
          << "\n";

  auto pipes_map = mesh_->links()->pipes();
  for (const auto& link : pipes_map) {
    auto pipe = link.second;
    auto node1 = pipe->nodes().first->id();
    auto node2 = pipe->nodes().second->id();
    auto pipe_prop = pipe->property();

    outfile << std::setprecision(12) << pipe_prop.name << "   " << node1
            << "   " << node2 << "   " << from_si(pipe_prop.length, "length")
            << "   " << from_si(pipe_prop.diameter, "diameter") << "   "
            << pipe_prop.roughness << "   " << pipe_prop.minor_loss_coeff
            << "   "
            << "Open"
            << "   "
            << "        ;"
            << "\n";
  }
}

void pipenetwork::IO_utils::Output::save_sim_result(
    const std::string& output_path) {
  create_new_folder(output_path);
  std::ofstream outnode(output_path + mesh_->name() + "_nodes.csv");
  std::ofstream outlink(output_path + mesh_->name() + "_links.csv");
  outnode << "node_name"
          << ","
          << "head"
          << ","
          << "demand"
          << "\n";
  outlink << "link_name"
          << ","
          << "flowrate"
          << "\n";

  // junctions
  auto junction_map = mesh_->nodes()->junctions();
  for (auto& index_junc : junction_map) {
    auto nid = index_junc.first;
    auto junction = index_junc.second;
    auto junc_prop = junction->property();
    outnode << std::setprecision(12) << junc_prop.name << "," << junction->head
            << "," << junction->demand << "\n";
  }
  // reservoirs
  auto res_map = mesh_->nodes()->reservoirs();
  for (auto& index_res : res_map) {
    auto nid = index_res.first;
    auto res = index_res.second;
    auto res_prop = res->property();
    outnode << std::setprecision(12) << res_prop.name << "," << res_prop.head
            << "," << res->discharge << "\n";
  }

  auto pipe_map = mesh_->links()->pipes();
  for (auto& index_pipe : pipe_map) {
    auto nid = index_pipe.first;
    auto pipe = index_pipe.second;
    auto pipe_prop = pipe->property();
    outlink << std::setprecision(12) << pipe_prop.name << "," << pipe->flowrate
            << "\n";
  }

  std::cout << "Simultation results saved to: " << output_path << std::endl;
}
