#include "catch.hpp"

#include "mesh.h"
#include "node.h"
#include "pipe.h"

// Check mesh class
TEST_CASE("Mesh is checked", "[Mesh]") {
  // Tolerance
  const double tolerance = 1.e-12;

  // Node indices
  std::array<const unsigned, 5> nodeids = {101, 102, 103, 104, 105};

  // Nodal coordinates
  const Eigen::Vector3d coords1(0.0, 0.0, 0.0);
  const Eigen::Vector3d coords2(2.0, 0.0, 0.0);
  const Eigen::Vector3d coords3(1.0, 1.0, 0.0);
  const Eigen::Vector3d coords4(0.0, 2.0, 0.0);
  const Eigen::Vector3d coords5(2.0, 2.0, 0.0);
  std::array<const Eigen::Vector3d, 5> coords = {coords1, coords2, coords3,
                                                 coords4, coords5};

  // Creat vector of nodal pointers
  std::vector<std::shared_ptr<Node>> nodes;
  for (unsigned i = 0; i < nodeids.size(); ++i)
    nodes.emplace_back(std::make_shared<Node>(nodeids.at(i), coords.at(i)));

  // Creat arrays of two node pointers with previous defined index and
  // coordinates of the nodes
  std::array<std::array<std::shared_ptr<Node>, 2>, 4> node_pair;
  node_pair.at(0) = {nodes.at(0), nodes.at(2)};
  node_pair.at(1) = {nodes.at(1), nodes.at(2)};
  node_pair.at(2) = {nodes.at(2), nodes.at(3)};
  node_pair.at(3) = {nodes.at(2), nodes.at(4)};

  // Pipe indices
  std::array<const unsigned, 4> pipeids = {201, 202, 203, 204};

  // Creat pipes based on previous created node pointers
  std::vector<std::shared_ptr<Pipe>> pipes;
  for (unsigned i = 0; i < pipeids.size(); ++i)
    pipes.emplace_back(std::make_unique<Pipe>(pipeids.at(i), node_pair.at(i)));

  // Mesh index
  const unsigned meshid = 301;

  // Creat a mesh
  auto mesh = std::make_unique<Mesh>(meshid);

  // Add nodes and pipes into the mesh
  for (unsigned i = 0; i < nodes.size(); ++i) mesh->add_node(nodes.at(i));
  for (unsigned i = 0; i < pipes.size(); ++i) mesh->add_pipe(pipes.at(i));

  // Check mesh id
  REQUIRE(mesh->id() == meshid);

  // Check redundant node
  mesh->redundant_node();

  // Check number of nodes
  REQUIRE(mesh->nnodes() == 5);

  // Check number of pipes
  REQUIRE(mesh->npipes() == 4);

  // Ckeck degree centrality of the nods
  REQUIRE(mesh->degree_centrality(0) == 1);
  REQUIRE(mesh->degree_centrality(2) == 4);

  // Check coordinates of nodes
  REQUIRE(mesh->nodal_coordinates().at(1)(0) == Approx(2.0).epsilon(tolerance));
  REQUIRE(mesh->nodal_coordinates().at(2)(0) == Approx(1.0).epsilon(tolerance));
}
