#include "catch.hpp"

#include "io.h"
#include "mesh.h"

// Check mesh class
TEST_CASE("Mesh is checked", "[Mesh]") {
  // Tolerance
  const double tolerance = 1.e-12;

  // create junctions
  std::vector<std::string> junction_names{"1", "2", "3", "4", "5"};
  std::vector<double> elevations{5, 4, 3, 2, 1};
  std::vector<double> demands{1, 2, 3, 4, 5};
  std::vector<double> leak_diameters{0, 0, 0, 0.2, 0.1};

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
  // 1<->2<->6; 4<->7

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
    REQUIRE(mesh_nodes.njunctions() == junc_props.size());
    // check reservoirs
    REQUIRE(mesh_nodes.nreservoirs() == res_props.size());
    // check nodes
    REQUIRE(mesh_nodes.nnodes() == (res_props.size() + junc_props.size()));

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
    REQUIRE(mesh_links.npipes() == pipe_props.size());
    // check pumps
    REQUIRE(mesh_links.npumps() == 0);
    // check links
    REQUIRE(mesh_links.nlinks() == pipe_props.size());
  }

  SECTION("Mesh") {
    std::string mesh_name = "test_mesh";
    auto mesh = pipenetwork::Mesh(mesh_name);
    mesh.create_nodes(junc_props, res_props);
    mesh.create_links(pipe_props, pump_props, valve_props);
    mesh.create_mesh_graph();
    // nodes
    auto mesh_nodes = mesh.nodes();
    REQUIRE(mesh_nodes->njunctions() == junc_props.size());

    // links
    auto mesh_links = mesh.links();
    // check pipes
    REQUIRE(mesh_links->npipes() == pipe_props.size());
    // check pumps
    REQUIRE(mesh_links->npumps() == 0);

    // check leak node ids
    auto leak_nids = mesh.leak_nids();
    REQUIRE(leak_nids.size() == 2);
    REQUIRE(leak_nids[0] == 3);
    REQUIRE(leak_nids[1] == 4);

    // graph
    auto graph = mesh.mesh_graph();
    // Adjacency matrix
    auto A = graph->adjacency_matrix();

    REQUIRE(A.coeff(0, 1) == 1);
    REQUIRE(A.coeff(1, 0) == 1);
    REQUIRE(A.coeff(1, 5) == 1);
    REQUIRE(A.coeff(5, 1) == 1);
    REQUIRE(A.coeff(3, 6) == 1);
    REQUIRE(A.coeff(6, 3) == 1);
    REQUIRE(A.coeff(4, 4) == 0);
    REQUIRE(A.coeff(2, 6) == 0);

    // Node to link map
    auto n2l_map = graph->node2link_map();
    REQUIRE(n2l_map[0].size() == 1);
    REQUIRE(n2l_map[1].size() == 2);
    REQUIRE(n2l_map[2].size() == 0);
    REQUIRE(n2l_map[3].size() == 1);
    REQUIRE(n2l_map[4].size() == 0);

    // Node degrees
    auto node_degrees = graph->ndegree();
    REQUIRE(node_degrees[0] == 1);
    REQUIRE(node_degrees[1] == 2);
    REQUIRE(node_degrees[2] == 0);

    // Bfs
    auto connectivity_mask1 = graph->bfs(1);
    REQUIRE(connectivity_mask1[0] == 1);
    REQUIRE(connectivity_mask1[1] == 1);
    REQUIRE(connectivity_mask1[5] == 1);
    REQUIRE(connectivity_mask1[3] == 0);
    REQUIRE(connectivity_mask1[6] == 0);
    auto connectivity_mask3 = graph->bfs(3);
    REQUIRE(connectivity_mask3[0] == 0);
    REQUIRE(connectivity_mask3[1] == 0);
    REQUIRE(connectivity_mask3[5] == 0);
    REQUIRE(connectivity_mask3[3] == 1);
    REQUIRE(connectivity_mask3[6] == 1);

    // isolated junctions
    auto iso_juncs = mesh.iso_nodes();
    REQUIRE(iso_juncs.size() == 2);
    REQUIRE(iso_juncs[0] == 2);
    REQUIRE(iso_juncs[1] == 4);

    // isolated links
    auto iso_links = mesh.iso_links();
    REQUIRE(iso_links.size() == 0);
  }

  SECTION("Mesh from IO") {
    auto IO = std::make_shared<pipenetwork::IO>();
    IO->read_inp("../test_files/test_net.inp");
    std::string mesh_name = "test_mesh_inp";
    auto mesh = pipenetwork::Mesh(mesh_name);
    mesh.create_nodes(IO->junction_properties(), IO->reservoir_properties());
    mesh.create_links(IO->pipe_properties(), IO->pump_properties(),
                      IO->valve_properties());
    mesh.create_mesh_graph();
    //    mesh.print_summary();

    REQUIRE(mesh.nodes()->njunctions() == 9);
    REQUIRE(mesh.nodes()->nreservoirs() == 2);
    REQUIRE(mesh.links()->npipes() == 8);
  }
}
