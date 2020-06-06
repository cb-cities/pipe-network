#include "catch.hpp"

#include "mesh.h"

// Check mesh class
TEST_CASE("Mesh is checked", "[Mesh]") {
  // Tolerance
  const double tolerance = 1.e-12;

  // create junctions
  std::vector<std::string> junction_names{"1", "2", "3", "4", "5"};
  std::vector<double> elevations{5, 4, 3, 2, 1};
  std::vector<double> demands{1, 2, 3, 4, 5};
  std::vector<double> leak_diameters{0.1, 0.4, 0.3, 0.2, 0.1};

  std::vector<pipenetwork::JunctionProp> junc_props;
  for (int i = 0; i < elevations.size(); ++i) {
    pipenetwork::JunctionProp junc_prop;
    junc_prop.name = junction_names[i];
    junc_prop.elevation = elevations[i];
    junc_prop.demand = demands[i];
    junc_prop.leak_diameter = leak_diameters[i];
    junc_props.emplace_back(junc_prop);
  }

  // Reservoirs
  std::vector<std::string> res_names{"6", "7"};
  std::vector<double> heads{99, 100};

  std::vector<pipenetwork::ReservoirProp> res_props;
  for (int i = 0; i < res_names.size(); ++i) {
    pipenetwork::ReservoirProp res_prop;
    res_prop.name = res_names[i];
    res_prop.head = heads[i];
    res_props.emplace_back(res_prop);
  }

  std::vector<std::string> pipe_names{"1", "2", "3"};
  std::vector<std::pair<std::string, std::string>> node_names{
      std::make_pair("1", "2"), std::make_pair("2", "6"),
      std::make_pair("4", "7")};
  const std::vector<double> length{100, 200, 300};
  const std::vector<double> diameter{3, 4, 5};
  const std::vector<double> roughness{0.2, .6, .9};

  std::vector<pipenetwork::PipeProp> pipe_props;
  for (int i = 0; i < pipe_names.size(); ++i) {
    pipenetwork::PipeProp pipe_prop;
    pipe_prop.name = pipe_names[i];
    pipe_prop.length = length[i];
    pipe_prop.diameter = diameter[i];
    pipe_prop.roughness = roughness[i];
    pipe_prop.node1_name = node_names[i].first;
    pipe_prop.node2_name = node_names[i].second;
    pipe_props.emplace_back(pipe_prop);
  }

  std::vector<pipenetwork::PumpProp> pump_props;
  std::vector<pipenetwork::ValveProp> valve_props;

  for (int i = 0; i < 3; ++i) {
    pipenetwork::PumpProp pump_prop;
    pump_prop.name = "pump-" + std::to_string(i);
    pump_props.emplace_back(pump_prop);
  }

  SECTION("MeshNodes") {
    auto mesh_nodes = pipenetwork::MeshNodes(junc_props, res_props);
    // check junctions
    REQUIRE(mesh_nodes.junctions.size() == junc_props.size());
    // check reservoirs
    REQUIRE(mesh_nodes.reservoirs.size() == res_props.size());

    // check get node
    auto junc1 = mesh_nodes.get_node("1");
    REQUIRE(junc1->id() == 0);

    auto res1 = mesh_nodes.get_node("6");
  }

  SECTION("MeshLinks") {
    auto mesh_nodes = pipenetwork::MeshNodes(junc_props, res_props);
    auto mesh_links =
        pipenetwork::MeshLinks(pipe_props, pump_props, valve_props, mesh_nodes);

    // check pipes
    REQUIRE(mesh_links.pipes.size() == pipe_props.size());
    // check pumps
    REQUIRE(mesh_links.pumps.size() == 0);
  }
  SECTION("Mesh") {
    std::string mesh_name = "test_mesh";
    auto mesh = pipenetwork::Mesh(mesh_name);
    mesh.create_nodes(junc_props, res_props);
    mesh.create_links(pipe_props, pump_props, valve_props);

    // nodes
    auto mesh_nodes = mesh.nodes();

    REQUIRE(mesh_nodes->junctions.size() == junc_props.size());

    // links
    auto mesh_links = mesh.links();
    // check pipes
    REQUIRE(mesh_links->pipes.size() == pipe_props.size());
    // check pumps
    REQUIRE(mesh_links->pumps.size() == 0);
  }

  //  std::vector<std::string> pipe_ids{"1", "2", "3"};
  //  std::vector<std::pair<std::string, std::string>> nodeids{
  //      std::make_pair("1", "2"), std::make_pair("2", "6"),
  //      std::make_pair("4", "7")};
  //  const std::vector<double> length{100, 200, 300};
  //  const std::vector<double> diameter{3, 4, 5};
  //  const std::vector<double> roughness{0.2, .6, .9};
  //  const std::vector<pipenetwork::Link_status> status{
  //      pipenetwork::OPEN, pipenetwork::OPEN, pipenetwork::OPEN};
  //
  //  std::vector<pipenetwork::Pipe_prop> pipe_props;
  //  for (int i = 0; i < pipe_ids.size(); ++i) {
  //    pipenetwork::Pipe_prop pipe_prop;
  //    pipe_prop.id = pipe_ids[i];
  //    pipe_prop.length = length[i];
  //    pipe_prop.diameter = diameter[i];
  //    pipe_prop.roughness = roughness[i];
  //    pipe_prop.node1_id = nodeids[i].first;
  //    pipe_prop.node2_id = nodeids[i].second;
  //    pipe_prop.status = status[i];
  //
  //    pipe_props.emplace_back(pipe_prop);
  //  }
  //
  //  mesh->create_pipes(pipe_props);
  //
  //  REQUIRE(mesh->njunctions() == 5);
  //  REQUIRE(mesh->nsources() == 2);
  //  REQUIRE(mesh->npipes() == 3);
  //  REQUIRE(mesh->nnodes() == 7);
  //
  //  auto node_map = mesh->nodes();
  //  auto links = mesh->links();
  //  auto junction_node = node_map.at("2");
  //  auto res_node = node_map.at("6");
  //  auto pipe = links[0];
  //
  //  REQUIRE(junction_node->nodal_info()["elevation"] == 4);
  //  REQUIRE(res_node->nodal_info()["head"] == 99);
  //  REQUIRE(pipe->link_info()["length"] == 100);

  //
  ////    mesh->print_summary ();
  //    double init_head = 10;
  //    mesh->iterate_over_nodes
  //    (std::bind(&pipenetwork::Node::update_sim_demand,
  //    std::placeholders::_1,init_head));
  //
}
