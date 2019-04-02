#include "io.h"

// Read input CSV file and store the data
bool pipenetwork::IO::read_network(const std::string& node_filename,
                                   const std::string& pipe_filename) {
  bool status = true;
  try {
    io::CSVReader<6> in_node(node_filename);
    io::CSVReader<7> in_pipe(pipe_filename);
    in_node.read_header(io::ignore_extra_column, "Node_id", "Node_coord1",
                        "Node_coord2", "Node_coord3", "Node_head",
                        "Node_discharge");
    in_pipe.read_header(io::ignore_extra_column, "Pipe_id", "Pipe_node1",
                        "Pipe_node2", "Pipe_diameter", "Pipe_roughness",
                        "Pipe_status", "Pipe_discharge");
    Index node_id, pipe_id, pipe_node1, pipe_node2;
    std::string node_head, node_discharge, pipe_discharge;
    long double node_coord1, node_coord2, node_coord3, pipe_diameter,
        pipe_roughness;
    int pipe_status;
    while (in_node.read_row(node_id, node_coord1, node_coord2, node_coord3,
                            node_head, node_discharge)) {
      Eigen::Vector3d coords(node_coord1, node_coord2, node_coord3);
      nodal_coords_.emplace_back(coords);

      if (node_head != "NA") {
        long double head = std::stod(node_head);
        initial_nodal_head_.emplace_back(std::make_pair(node_id, head));
      }
      if (node_discharge != "NA") {
        long double ndischarge = std::stod(node_discharge);
        initial_nodal_discharge_.emplace_back(
            std::make_pair(node_id, ndischarge));
      }
    }

    while (in_pipe.read_row(pipe_id, pipe_node1, pipe_node2, pipe_diameter,
                            pipe_roughness, pipe_status, pipe_discharge)) {
      node_pairs_.emplace_back(std::make_pair(pipe_node1, pipe_node2));
      diameters_.emplace_back(pipe_diameter);
      roughness_.emplace_back(pipe_roughness);
      if (pipe_status == 0)
        pipe_status_.emplace_back(false);
      else
        pipe_status_.emplace_back(true);

      if (pipe_discharge != "NA") {
        long double pdischarge = std::stod(pipe_discharge);
        initial_pipe_discharge_.emplace_back(
            std::make_pair(pipe_id, pdischarge));
      }
    }

    std::cout << "Read from: " << node_filename << " and " << pipe_filename
              << std::endl
              << "The network contains " << nodal_coords_.size()
              << " nodes and " << node_pairs_.size() << "pipes." << std::endl
              << "With " << initial_nodal_head_.size()
              << " known initial nodal heads, "
              << initial_nodal_discharge_.size()
              << " known initial nodal discharges and "
              << initial_pipe_discharge_.size()
              << " known initial pipe discharges." << std::endl;
  } catch (std::exception& exception) {
    std::cout << "Read pipe network file: " << exception.what() << std::endl;
    status = false;
  }

  return status;
}
