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
  auto mesh = std::make_shared<Mesh>(meshid);

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
  std::vector<std::pair<Index, Index>> node_pairs;
  node_pairs.emplace_back(std::make_pair(0, 2));
  node_pairs.emplace_back(std::make_pair(1, 2));
  node_pairs.emplace_back(std::make_pair(2, 3));
  node_pairs.emplace_back(std::make_pair(2, 4));

  // Input vector of pipe diameter and status
  std::vector<double> diameter{1, 1, 1, 1};
  std::vector<double> roughness{0.01, 0.01, 0.01, 0.01};
  std::vector<bool> status{true, true, true, true};

  // Check mesh id
  REQUIRE(mesh->id() == meshid);

  // Check the case in which exception is not thrown
  SECTION("Check the case in which exception is not thrown") {

    // Create pipes based on pipe indices and previous created node pointers in
    // the mesh
    bool all_pipe_created =
        mesh->create_pipes(node_pairs, diameter, roughness, status);

    // Check number of nodes before remove
    REQUIRE(mesh->nnodes() == 7);

    // Remove isolated nodes from mesh
    mesh->remove_unconnected_nodes();

    // Check number of nodes after remove
    REQUIRE(mesh->nnodes() == 5);
    // Check number of pipes
    REQUIRE(mesh->npipes() == 4);
    // Check whether exception is thrown when creating pipe
    REQUIRE(all_pipe_created == true);
  }

  // Check the case in which exception is thrown
  SECTION("Check the case in which exception is thrown") {

    // Make pairs of nodes which do not exist to create pipe
    node_pairs.emplace_back(std::make_pair(2, 8));
    node_pairs.emplace_back(std::make_pair(8, 9));
    std::vector<double> diameter{1, 1, 1, 1, 1, 1};
    std::vector<double> roughness{0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
    std::vector<bool> status{true, true, true, true, true, true};
    // Create pipes based on pipe indices and previous created node pointers in
    // the mesh
    bool all_pipe_created =
        mesh->create_pipes(node_pairs, diameter, roughness, status);

    // Check number of nodes before remove
    REQUIRE(mesh->nnodes() == 7);

    // Remove isolated nodes from mesh
    mesh->remove_unconnected_nodes();

    // Check number of nodes after remove
    REQUIRE(mesh->nnodes() == 5);
    // Check number of pipes
    REQUIRE(mesh->npipes() == 4);
    // Check whether exception is thrown when creating pipe
    REQUIRE(all_pipe_created == false);
  }
}
