#include "catch.hpp"

#include "junction.h"
#include "mesh.h"
#include "pipe.h"
#include "reservoir.h"

// Check mesh class
TEST_CASE("Mesh is checked", "[Mesh]") {
  // Tolerance
  const double tolerance = 1.e-12;

  // Mesh index
  std::string meshid = "Testing Net";

  // Creat a mesh
  auto mesh = std::make_unique<pipenetwork::Mesh>(meshid);

  // Junctions
  pipenetwork::Reservoir_prop res_1;
  res_1.id = "321";
  res_1.head = 10;

  pipenetwork::Junction_prop junction_1;
  junction_1.id = "123";
  junction_1.elevation = 10;
  junction_1.demand = 20;
  junction_1.leak_diameter = 5;

  std::vector<std::string> junction_ids{"1", "2", "3", "4", "5"};
  std::vector<double> elevations{5, 4, 3, 2, 1};
  std::vector<double> demands{1, 2, 3, 4, 5};
  std::vector<double> leak_diameters{0.1, 0.4, 0.3, 0.2, 0.1};

  std::vector<pipenetwork::Junction_prop> junc_props;
  for (int i = 0; i < elevations.size(); ++i) {
    pipenetwork::Junction_prop junc_prop;
    junc_prop.id = junction_ids[i];
    junc_prop.elevation = elevations[i];
    junc_prop.demand = demands[i];
    junc_prop.leak_diameter = leak_diameters[i];

    junc_props.emplace_back(junc_prop);
  }

  mesh->create_junctions(junc_props);

  // Reservoirs
  std::vector<std::string> res_ids{"6", "7"};
  std::vector<double> heads{99, 100};

  std::vector<pipenetwork::Reservoir_prop> res_props;
  for (int i = 0; i < res_ids.size(); ++i) {
    pipenetwork::Reservoir_prop res_prop;
    res_prop.id = res_ids[i];
    res_prop.head = heads[i];

    res_props.emplace_back(res_prop);
  }
  mesh->create_reservoirs(res_props);

  std::vector<std::string> pipe_ids{"1", "2", "3"};
  std::vector<std::pair<std::string, std::string>> nodeids{
      std::make_pair("1", "2"), std::make_pair("2", "6"),
      std::make_pair("4", "7")};
  const std::vector<double> length{100, 200, 300};
  const std::vector<double> diameter{3, 4, 5};
  const std::vector<double> roughness{0.2, .6, .9};
  const std::vector<pipenetwork::Link_status> status{pipenetwork::OPEN, pipenetwork::OPEN, pipenetwork::OPEN};

  std::vector<pipenetwork::Pipe_prop> pipe_props;
  for (int i = 0; i < pipe_ids.size(); ++i) {
    pipenetwork::Pipe_prop pipe_prop;
    pipe_prop.id = pipe_ids[i];
    pipe_prop.length = length[i];
    pipe_prop.diameter = diameter[i];
    pipe_prop.roughness = roughness[i];
    pipe_prop.node1_id = nodeids[i].first;
    pipe_prop.node2_id = nodeids[i].second;
    pipe_prop.status = status[i];

    pipe_props.emplace_back(pipe_prop);
  }

  mesh->create_pipes(pipe_props);
  //  mesh->print_summary();

  //
  ////    mesh->print_summary ();
  //    double init_head = 10;
  //    mesh->iterate_over_nodes
  //    (std::bind(&pipenetwork::Node::update_sim_demand,
  //    std::placeholders::_1,init_head));
  //
}
