#include "catch.hpp"

#include "mesh.h"
#include "node.h"
#include "pipe.h"
#include "settings.h"

// Check mesh class
TEST_CASE("Mesh is checked", "[Mesh]") {
  // Tolerance
  const double tolerance = 1.e-12;

  // Mesh index
  const unsigned meshid = 101;

  // Creat a mesh
  auto mesh = std::make_unique<Mesh>(meshid);

  // Nodal coordinates
  const Eigen::Vector3d coords1(0.0, 0.0, 0.0);
  const Eigen::Vector3d coords2(2.0, 0.0, 0.0);
  const Eigen::Vector3d coords3(1.0, 1.0, 0.0);
  const Eigen::Vector3d coords4(0.0, 2.0, 0.0);
  const Eigen::Vector3d coords5(2.0, 2.0, 0.0);
  const Eigen::Vector3d coords6(3.0, 3.0, 0.0);
  const Eigen::Vector3d coords7(3.0, 4.0, 0.0);
  std::vector<Eigen::Vector3d> coords = {coords1, coords2, coords3, coords4,
                                         coords5, coords6, coords7};

  // Create nodal pointers based on nodal coordinates in the mesh
  mesh->create_nodes(coords);

  // Make pairs of nodes to create pipe
  std::vector<std::pair<Index, Index>> nodepair;
  nodepair.emplace_back(std::make_pair(0, 2));
  nodepair.emplace_back(std::make_pair(1, 2));
  nodepair.emplace_back(std::make_pair(2, 3));
  nodepair.emplace_back(std::make_pair(2, 4));

  // Create pipes based on pipe indices and previous created node pointers in
  // the mesh
  mesh->create_pipes(nodepair);

  // Remove isolated nodes from mesh and record them
  std::vector<std::shared_ptr<pipenetwork::Node>> isolated_nodes =
      mesh->isolated_nodes();

  // Check mesh id
  REQUIRE(mesh->id() == meshid);

  // Check isolated node
  // REQUIRE(mesh->isolated_node() == false);

  // Check number of nodes
  REQUIRE(mesh->nnodes() == 5);

  // Check number of pipes
  REQUIRE(mesh->npipes() == 4);

  // Check isolated nodes
  REQUIRE(isolated_nodes.size() == 2);
  REQUIRE(isolated_nodes.at(0)->coordinates()(1) ==
          Approx(3.0).epsilon(tolerance));
  REQUIRE(isolated_nodes.at(1)->coordinates()(1) ==
          Approx(4.0).epsilon(tolerance));
}
