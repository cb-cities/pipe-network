#include "catch.hpp"

#include "mesh.h"
#include "node.h"
#include "pipe.h"

// Check mesh class
TEST_CASE("Mesh is checked", "[Mesh]") {
  // Tolerance
  const double tolerance = 1.e-12;

  // Mesh index
  const unsigned meshid = 101;

  // Creat a mesh
  auto mesh = std::make_unique<Mesh>(meshid);

  // Node indices
  std::vector<unsigned> nodeids = {201, 202, 203, 204, 205};

  // Nodal coordinates
  const Eigen::Vector3d coords1(0.0, 0.0, 0.0);
  const Eigen::Vector3d coords2(2.0, 0.0, 0.0);
  const Eigen::Vector3d coords3(1.0, 1.0, 0.0);
  const Eigen::Vector3d coords4(0.0, 2.0, 0.0);
  const Eigen::Vector3d coords5(2.0, 2.0, 0.0);
  std::vector<Eigen::Vector3d> coords = {coords1, coords2, coords3, coords4,
                                         coords5};

  // Create nodal pointers based on nodal indices and coordinates in the mesh
  mesh->create_node(nodeids, coords);

  // Pipe indices
  std::vector<unsigned> pipeids = {301, 302, 303, 304};

  // Make pairs of nodes to create pipe
  std::vector<std::pair<unsigned, unsigned>> nodepair;
  nodepair.emplace_back(std::make_pair(nodeids.at(0), nodeids.at(2)));
  nodepair.emplace_back(std::make_pair(nodeids.at(1), nodeids.at(2)));
  nodepair.emplace_back(std::make_pair(nodeids.at(2), nodeids.at(3)));
  nodepair.emplace_back(std::make_pair(nodeids.at(2), nodeids.at(4)));

  // Create pipes based on pipe indices and previous created node pointers in
  // the mesh
  mesh->create_pipe(pipeids, nodepair);

  // Check mesh id
  REQUIRE(mesh->id() == meshid);

  // Check isolated node
  REQUIRE(mesh->isolated_node() == false);

  // Check number of nodes
  REQUIRE(mesh->nnodes() == 5);

  // Check number of pipes
  REQUIRE(mesh->npipes() == 4);

  // Check coordinates of nodes
  REQUIRE(mesh->nodal_coordinates().at(1)(0) == Approx(2.0).epsilon(tolerance));
  REQUIRE(mesh->nodal_coordinates().at(2)(0) == Approx(1.0).epsilon(tolerance));
}
