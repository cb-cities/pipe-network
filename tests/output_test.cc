#include "catch.hpp"

#include "input.h"
#include "mesh.h"
#include "output.h"
// Check output class
TEST_CASE("Output is checked", "[Output]") {
  double tolerance = 1e-6;
  int syn_size = 300;
  auto Input = std::make_shared<pipenetwork::Input>(syn_size);

  // Mesh ID
  std::string meshid = "syn_net_" + std::to_string(syn_size);
  // Creat a mesh
  auto mesh = std::make_shared<pipenetwork::Mesh>(meshid);
  mesh->create_mesh_from_inp(Input);

  auto Output = std::make_shared<pipenetwork::Output>(
      mesh, "../benchmarks/write_test_large.inp");
}