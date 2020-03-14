#include "output.h"

void pipenetwork::Output::write_junctions() {
  // write row info
  outFile_ << "[JUNCTIONS]"
           << "\n";
  outFile_ << ";ID    "
           << "Elev   "
           << "Demand     "
           << "Pattern     "
           << "\n";
  for (const auto& node : mesh_->nodes()) {
    // get node information, determine possible leak nodes and assemble demands
    // heads vector
    auto info = node.second->nodal_info();
    auto junction_id = node.second->id();
    switch (static_cast<int>(info["type"])) {
      case JUNCTION:
        outFile_ << std::setprecision(12) << junction_id << "   "
                 << from_si(info["elevation"], "elevation") << "   "
                 << from_si(info["demand"], "demand") << "   "
                 << "        ;"
                 << "\n";
        break;
    }
  }
}

void pipenetwork::Output::write_reservoirs() {
  // write row info
  outFile_ << "\n";
  outFile_ << "[RESERVOIRS]"
           << "\n";
  outFile_ << ";ID    "
           << "Head   "
           << "Pattern     "
           << "\n";
  for (const auto& node : mesh_->nodes()) {
    // get node information, determine possible leak nodes and assemble demands
    // heads vector
    auto info = node.second->nodal_info();
    auto junction_id = node.second->id();
    switch (static_cast<int>(info["type"])) {
      case RESERVOIR:
        outFile_ << std::setprecision(12) << junction_id << "   "
                 << from_si(info["head"], "head") << "   "
                 << "        ;"
                 << "\n";
        break;
    }
  }
}

void pipenetwork::Output::write_pipes() {
  // write row info
  outFile_ << "\n";
  outFile_ << "[PIPES]"
           << "\n";
  outFile_ << ";ID    "
           << "Node1   "
           << "Node2   "
           << "Length   "
           << "Diameter   "
           << "Roughness   "
           << "MinorLoss   "
           << "Status   "
           << "\n";
  std::string status;
  for (const auto& link : mesh_->links()) {
    // get node information, determine possible leak nodes and assemble demands
    // heads vector
    auto info = link->link_info();
    auto link_id = link->id();
    auto node1 = link->nodes().first->id();
    auto node2 = link->nodes().second->id();
    if (link->link_status() == OPEN) {
      status = "Open";
    } else {
      status = "Closed";
    }

    switch (static_cast<int>(info["type"])) {
      case PIPE:
        outFile_ << std::setprecision(12) << link_id << "   " << node1 << "   "
                 << node2 << "   " << from_si(info["length"], "length") << "   "
                 << from_si(info["diameter"], "diameter") << "   "
                 << info["roughness"] << "   " << info["minor_loss"] << "   "
                 << status << "   "
                 << "        ;"
                 << "\n";
        break;
    }
  }
}

double pipenetwork::from_si(double val, const std::string& mode) {
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
