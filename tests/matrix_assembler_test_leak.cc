#include "catch.hpp"
#include <cmath>

#include "matrix_assembler.h"
#include "mesh.h"
#include "node.h"
#include "pipe.h"
#include "settings.h"

// Check matrix_assembler class
TEST_CASE("MatrixAssembler is checked", "[MatrixAssembler]") {
  // Tolerance
  const double tolerance = 1.e-12;

  // Mesh index
  const unsigned meshid = 102;

  // Creat a mesh
  auto mesh = std::make_shared<pipenetwork::Mesh>(meshid);

  // Nodal coordinates
  // using example netwrok in Todini(2013)
  const Eigen::Vector3d coords1(1.0, 3.0, 0.0);
  const Eigen::Vector3d coords2(0.0, 2.0, 0.0);
  const Eigen::Vector3d coords3(2.0, 2.0, 0.0);
  const Eigen::Vector3d coords4(0.0, 0.0, 0.0);
  const Eigen::Vector3d coords5(2.0, 0.0, 0.0);
  std::vector<Eigen::Vector3d> coords = {coords1, coords2, coords3, coords4,
                                         coords5};

  // Create nodal pointers based on nodal coordinates in the mesh
  mesh->create_nodes(coords);

  // Make pairs of nodes to create pipe
  std::vector<std::pair<Index, Index>> node_pairs;
  node_pairs.emplace_back(std::make_pair(0, 1));
  node_pairs.emplace_back(std::make_pair(0, 2));
  node_pairs.emplace_back(std::make_pair(1, 2));
  node_pairs.emplace_back(std::make_pair(1, 3));
  node_pairs.emplace_back(std::make_pair(1, 4));
  node_pairs.emplace_back(std::make_pair(2, 4));
  node_pairs.emplace_back(std::make_pair(3, 4));

  // Input vector of pipe diameter, roughness and status
  std::vector<double> diameter{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};
  std::vector<double> roughness{
      147683.36372579, 4102.7779645159, 1238.1607371737, 1598.5610759603,
      430.5609768094,  137.58891975265, 133.37351955511};
  std::vector<bool> status{true, true, true, true, true, true, true};

  // Create pipes based on pipe indices and previous created node pointers in
  // the mesh
  bool all_pipe_created =
      mesh->create_pipes(node_pairs, diameter, roughness, status);

  // Assign initial nodal head and discharge
  double default_init_node_head = 0.0;
  double default_init_node_discharge = 10;
  double default_init_pipe_discharge = 0.001;
  double head1 = 100.0;
  double head2 = 99.0;
  double discharge1 = 100.0;
  double discharge2 = -10.0;
  double leak_dia1 = 0.1;
  double leak_dia2 = 0.2;

  // Initialize pipe discharge
  mesh->initialize_pipe_discharge(default_init_pipe_discharge);

  std::vector<std::pair<Index, double>> node_elevation;
  node_elevation.emplace_back(std::make_pair(0, head1));
  node_elevation.emplace_back(std::make_pair(1, head2));
  std::vector<std::pair<Index, double>> node_discharge;
  node_discharge.emplace_back(std::make_pair(0, discharge1));
  node_discharge.emplace_back(std::make_pair(1, discharge2));
  std::vector<std::pair<Index, double>> node_leak_dia;
  node_leak_dia.emplace_back(std::make_pair(0, leak_dia1));
  node_leak_dia.emplace_back(std::make_pair(1, leak_dia2));

  mesh->assign_node_elevation(node_elevation);
  mesh->assign_node_demand(node_discharge);
  mesh->assign_node_leak(node_leak_dia);

  SECTION("DD MODE") {
    // Initialize matrix assembler and obtain global index to nodes and pipes
    bool pdd_mode = false;
    auto assembler = std::make_shared<pipenetwork::MatrixAssembler>(pdd_mode);
    assembler->global_nodal_pipe_indices(mesh);

    // Check the number of nodes and pipes in the network
    unsigned nnodes = assembler->nnodes();
    unsigned npipes = assembler->npipes();
    REQUIRE(nnodes == 5);
    REQUIRE(npipes == 7);
    // Check initialized variable (head and discharge) vector
    SECTION("Check initialized variable vector") {
      // Initialize variable vector
      assembler->assemble_variable_vector();
      std::shared_ptr<Eigen::VectorXd> variable_vec = assembler->variable_vec();
      //      std::cout << (*variable_vec) << std::endl;
      // Check nodal head
      REQUIRE(variable_vec->coeff(0) == Approx(head1).epsilon(tolerance));
      REQUIRE(variable_vec->coeff(1) == Approx(head2).epsilon(tolerance));
      for (int i = 2; i <= 4; i++)
        REQUIRE(variable_vec->coeff(i) ==
                Approx(default_init_node_head).epsilon(tolerance));
      // Check nodal discharge
      REQUIRE(variable_vec->coeff(5) == Approx(discharge1).epsilon(tolerance));
      REQUIRE(variable_vec->coeff(6) == Approx(discharge2).epsilon(tolerance));
      for (int i = 7; i <= 9; i++)
        REQUIRE(variable_vec->coeff(i) ==
                Approx(default_init_node_discharge).epsilon(tolerance));
      // Check pipe discharge
      for (int i = 10; i <= 16; i++)
        REQUIRE(variable_vec->coeff(i) ==
                Approx(default_init_pipe_discharge).epsilon(tolerance));
    }

    // Check initialized residual vector
    SECTION("Check initialized residual vector") {
      // Initialize residual vector
      assembler->assemble_residual_vector();
      std::shared_ptr<Eigen::VectorXd> residual_vec = assembler->residual_vec();
      //                  std::cout << (*residual_vec) << std::endl;
      // Check nodal balance residual
      REQUIRE(residual_vec->coeff(0) == Approx(100.002).epsilon(tolerance));
      REQUIRE(residual_vec->coeff(1) == Approx(-9.998).epsilon(tolerance));
      REQUIRE(residual_vec->coeff(2) == Approx(9.999).epsilon(tolerance));
      REQUIRE(residual_vec->coeff(3) == Approx(10).epsilon(tolerance));
      REQUIRE(residual_vec->coeff(4) == Approx(9.997).epsilon(tolerance));
      // Check headloss residual
      REQUIRE(residual_vec->coeff(5) ==
              Approx(-0.99999999916924).epsilon(tolerance));
      REQUIRE(residual_vec->coeff(6) ==
              Approx(-99.999999978347).epsilon(tolerance));
      REQUIRE(residual_vec->coeff(7) ==
              Approx(-98.999999960916).epsilon(tolerance));
      REQUIRE(residual_vec->coeff(8) ==
              Approx(-98.999999994002).epsilon(tolerance));
      REQUIRE(residual_vec->coeff(9) ==
              Approx(-98.99999996752).epsilon(tolerance));
      REQUIRE(residual_vec->coeff(10) ==
              Approx(0.0000000781682).epsilon(tolerance));
      REQUIRE(residual_vec->coeff(11) ==
              Approx(0.0000000390841).epsilon(tolerance));
    }
    //         Check initialized Jacobian matrix
    SECTION("Check initialized Jacobian matrix") {

      // Assemble Jacobian matrix
      assembler->assemble_jacobian();
      std::shared_ptr<Eigen::SparseMatrix<double>> jac = assembler->jac();
      //              std::cout << (*jac) << std::endl;

      // Check jacB (nodal discharge in nodal balance equation)
      for (int row = 0; row < nnodes; row++) {
        for (int col = nnodes; col < nnodes + nnodes; col++) {
          if ((col - nnodes) == row) {
            REQUIRE(jac->coeff(row, col) == -1);
          } else {
            REQUIRE(jac->coeff(row, col) == 0);
          }
        }
      }

      // Check jacC (pipe discharge in nodal balance equation)
      REQUIRE(jac->coeff(0, 10) == -1);
      REQUIRE(jac->coeff(0, 11) == -1);
      REQUIRE(jac->coeff(1, 10) == 1);
      REQUIRE(jac->coeff(1, 12) == -1);
      REQUIRE(jac->coeff(1, 13) == -1);
      REQUIRE(jac->coeff(1, 14) == -1);
      REQUIRE(jac->coeff(2, 11) == 1);
      REQUIRE(jac->coeff(2, 12) == 1);
      REQUIRE(jac->coeff(2, 15) == -1);
      REQUIRE(jac->coeff(3, 13) == 1);
      REQUIRE(jac->coeff(3, 16) == -1);
      REQUIRE(jac->coeff(4, 14) == 1);
      REQUIRE(jac->coeff(4, 15) == 1);
      REQUIRE(jac->coeff(4, 16) == 1);

      // Check jacD (nodal head in headloss equation)
      REQUIRE(jac->coeff(5, 0) == 1);
      REQUIRE(jac->coeff(5, 1) == -1);
      REQUIRE(jac->coeff(6, 0) == 1);
      REQUIRE(jac->coeff(6, 2) == -1);
      REQUIRE(jac->coeff(7, 1) == 1);
      REQUIRE(jac->coeff(7, 2) == -1);
      REQUIRE(jac->coeff(8, 1) == 1);
      REQUIRE(jac->coeff(8, 3) == -1);
      REQUIRE(jac->coeff(9, 1) == 1);
      REQUIRE(jac->coeff(9, 4) == -1);
      REQUIRE(jac->coeff(10, 2) == 1);
      REQUIRE(jac->coeff(10, 4) == -1);
      REQUIRE(jac->coeff(11, 3) == 1);
      REQUIRE(jac->coeff(11, 4) == -1);

      // Check jacF (pipe discharge in headloss equation)
      REQUIRE(jac->coeff(5, 10) ==
              Approx(-0.000001538573761089).epsilon(tolerance));
      REQUIRE(jac->coeff(6, 11) ==
              Approx(-0.00004010175932682).epsilon(tolerance));
      REQUIRE(jac->coeff(7, 12) ==
              Approx(-0.00007238373389406).epsilon(tolerance));
      REQUIRE(jac->coeff(8, 13) ==
              Approx(-0.00001110851163928).epsilon(tolerance));
      REQUIRE(jac->coeff(9, 14) ==
              Approx(-0.00006015263899023).epsilon(tolerance));
      REQUIRE(jac->coeff(10, 15) ==
              Approx(-0.0001447674677881).epsilon(tolerance));
      REQUIRE(jac->coeff(11, 16) ==
              Approx(-0.00007238373389406).epsilon(tolerance));
      REQUIRE(jac->coeff(13, 1) == Approx(1).epsilon(tolerance));
    }
  }

  SECTION("PDD MODE") {
    // Initialize matrix assembler and obtain global index to nodes and pipes
    bool pdd_mode = true;
    auto assembler = std::make_shared<pipenetwork::MatrixAssembler>(pdd_mode);

    std::vector<std::pair<Index, double>> node_head;
    node_head.emplace_back(std::make_pair(0, 110));
    node_head.emplace_back(std::make_pair(1, 100));

    mesh->assign_node_head(node_head);
    assembler->global_nodal_pipe_indices(mesh);

    // Check the number of nodes and pipes in the network
    unsigned nnodes = assembler->nnodes();
    unsigned npipes = assembler->npipes();
    assembler->assemble_variable_vector();
    std::shared_ptr<Eigen::VectorXd> variable_vec = assembler->variable_vec();
    //    std::cout << (*variable_vec)<< std::endl;
    // Check initialized residual vector
    SECTION("Check initialized residual vector") {
      // Initialize residual vector
      assembler->assemble_residual_vector();
      std::shared_ptr<Eigen::VectorXd> residual_vec = assembler->residual_vec();
      //      std::cout << (*variable_vec) << std::endl;
      //      std::cout << (*residual_vec) << std::endl;

      // Check pressure demand residual
      double res =
          -1 * (100 - 100 * std::pow((110.0 - 100 - 0) / (20 - 0), 0.5));
      REQUIRE(residual_vec->coeff(12) == Approx(res).epsilon(tolerance));
      // Check pressure leak residual
      double leak1 = 0.75 * 3.1415926 * std::pow(leak_dia1 / 2, 2) *
                     std::pow(2 * 9.81 * 10, 0.5);
      double leak2 = 0.75 * 3.1415926 * std::pow(leak_dia2 / 2, 2) *
                     std::pow(2 * 9.81 * 1, 0.5);
      REQUIRE(residual_vec->coeff(17) == Approx(leak1).epsilon(tolerance));
      REQUIRE(residual_vec->coeff(18) == Approx(leak2).epsilon(tolerance));
    }
    //         Check initialized Jacobian matrix
    SECTION("Check initialized Jacobian matrix") {

      // Assemble Jacobian matrix
      assembler->assemble_jacobian();
      std::shared_ptr<Eigen::SparseMatrix<double>> jac = assembler->jac();
      //      std::cout << (*jac) << std::endl;

      // Check jacG (demand equation)
      double res = -1 * (0.5 * 100 * std::pow((10.0) / (20 - 0), -0.5));
      REQUIRE(jac->coeff(12, 0) == Approx(res).epsilon(tolerance));
      // check jacJ (leak equation)
      double leak_jac_1 = 0.5 * 0.75 * 3.1415926 * std::pow(leak_dia1 / 2, 2) *
                          std::pow(2 * 9.81, 0.5) * std::pow(10, -0.5);
      double leak_jac_2 = 0.5 * 0.75 * 3.1415926 * std::pow(leak_dia2 / 2, 2) *
                          std::pow(2 * 9.81, 0.5) * std::pow(1, -0.5);
      REQUIRE(jac->coeff(17, 0) == Approx(leak_jac_1).epsilon(tolerance));
      REQUIRE(jac->coeff(18, 1) == Approx(leak_jac_2).epsilon(tolerance));
    }
  }
}