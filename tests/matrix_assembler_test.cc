#include "catch.hpp"
#include <cmath>

#include "matrix_assembler.h"
// Check matrix_assembler class
TEST_CASE("MatrixAssembler is checked", "[MatrixAssembler]") {

  // Tolerance
  const double tolerance = 1.e-6;

  // Mesh index
  const unsigned meshid = 101;

  // Creat a mesh
  auto mesh = std::make_shared<pipenetwork::Mesh>(meshid);

  std::vector<Index> junction_ids{10, 11, 12};
  std::vector<double> elevations{216.408, 216.408, 213.36};
  std::vector<double> demands{0, 9.464e-03, 9.464e-03};
  std::vector<double> leak_diameters{0, 0.1, 0};

  mesh->create_junctions(junction_ids, elevations, demands, leak_diameters);

  std::vector<Index> res_ids{13};
  std::vector<double> heads{3.048e+02};

  mesh->create_reservoirs(res_ids, heads);

  std::vector<Index> pipe_ids{9, 10, 11, 12};
  std::vector<std::pair<Index, Index>> nodeids{
      std::make_pair(13, 10), std::make_pair(10, 11), std::make_pair(11, 12),
      std::make_pair(10, 12)};
  const std::vector<double> length{3209.5440000000003, 3209.5440000000003,
                                   1609.344, 1609.344};
  const std::vector<double> diameter{0.5588, 0.4572, 0.35559999999999997,
                                     0.254};
  const std::vector<double> roughness{100, 100, 100, 100};
  const std::vector<Pipe_status> status{OPEN, OPEN, OPEN, OPEN};

  mesh->create_pipes(pipe_ids, nodeids, length, diameter, roughness, status);

  //    mesh->print_summary ();

  double init_discharge = 1e-3;

  mesh->iterate_over_links(std::bind(&pipenetwork::Link::update_sim_discharge,
                                     std::placeholders::_1,
                                     init_discharge));  // initialze discharge
  SECTION("DD MODE") {
    bool pdd_mode = false;
    auto assembler =
        std::make_shared<pipenetwork::MatrixAssembler>(mesh, pdd_mode);

    SECTION("vairable test") {
      auto var_vec = (assembler->variable_vector());
      //      std::cout << (*var_vec) << std::endl;

      REQUIRE(var_vec->coeff(0) == Approx(216.408).epsilon(tolerance));
      REQUIRE(var_vec->coeff(5) == Approx(9.464e-03).epsilon(tolerance));
      REQUIRE(var_vec->coeff(10) == Approx(init_discharge).epsilon(tolerance));
      REQUIRE(var_vec->coeff(12) == Approx(0).epsilon(tolerance));
    }
    SECTION("residual test") {
      assembler->assemble_residual();
      auto res_vec = assembler->residual_vector();
//      std::cout << (*res_vec) << std::endl;
      REQUIRE(res_vec->coeff(1) == Approx(-9.464e-03).epsilon(tolerance));
      REQUIRE(res_vec->coeff(5) == Approx(0).epsilon(tolerance));
      REQUIRE(res_vec->coeff(8) == Approx(-88.3916796736).epsilon(tolerance));
      REQUIRE(res_vec->coeff(12) == Approx(0).epsilon(tolerance));
    }
    SECTION("jacobian test: static") {
      auto jac_matrix = assembler->jac_matrix();
      //      std::cout << (*jac_matrix) << std::endl;
      REQUIRE(jac_matrix->coeff(0, 4) == -1);
      REQUIRE(jac_matrix->coeff(1, 5) == -1);
      REQUIRE(jac_matrix->coeff(2, 6) == -1);
      REQUIRE(jac_matrix->coeff(3, 7) == -1);
      REQUIRE(jac_matrix->coeff(0, 9) == -1);
      REQUIRE(jac_matrix->coeff(0, 11) == -1);
      REQUIRE(jac_matrix->coeff(0, 8) == 1);
      REQUIRE(jac_matrix->coeff(1, 10) == -1);
      REQUIRE(jac_matrix->coeff(1, 9) == 1);
      REQUIRE(jac_matrix->coeff(2, 10) == 1);
      REQUIRE(jac_matrix->coeff(2, 11) == 1);
      REQUIRE(jac_matrix->coeff(3, 8) == -1);
      REQUIRE(jac_matrix->coeff(1, 12) == -1);
      REQUIRE(jac_matrix->coeff(4, 0) == 0);
      REQUIRE(jac_matrix->coeff(5, 1) == 0);
      REQUIRE(jac_matrix->coeff(6, 2) == 0);
      REQUIRE(jac_matrix->coeff(7, 3) == 1);
      REQUIRE(jac_matrix->coeff(4, 4) == 1);
      REQUIRE(jac_matrix->coeff(5, 5) == 1);
      REQUIRE(jac_matrix->coeff(6, 6) == 1);
      REQUIRE(jac_matrix->coeff(7, 7) == 0);
      REQUIRE(jac_matrix->coeff(8, 3) == -1);
      REQUIRE(jac_matrix->coeff(9, 0) == -1);
      REQUIRE(jac_matrix->coeff(10, 1) == -1);
      REQUIRE(jac_matrix->coeff(11, 0) == -1);
      REQUIRE(jac_matrix->coeff(8, 0) == 1);
      REQUIRE(jac_matrix->coeff(9, 1) == 1);
      REQUIRE(jac_matrix->coeff(10, 2) == 1);
      REQUIRE(jac_matrix->coeff(11, 2) == 1);
      REQUIRE(jac_matrix->coeff(8, 8) == 1);
      REQUIRE(jac_matrix->coeff(9, 9) == 1);
      REQUIRE(jac_matrix->coeff(10, 10) == 1);
      REQUIRE(jac_matrix->coeff(11, 11) == 1);
      REQUIRE(jac_matrix->coeff(12, 1) == 0);
      REQUIRE(jac_matrix->coeff(12, 12) == 1);
    }
    SECTION("jacobian test: dynamic, case 1 (discharge larger than HWQ2)") {
      auto jac_matrix = assembler->jac_matrix();
      assembler->update_jacobian();
//      std::cout << (*jac_matrix) << std::endl;
      REQUIRE(jac_matrix->coeff(0, 4) == -1);
      REQUIRE(jac_matrix->coeff(1, 5) == -1);
      REQUIRE(jac_matrix->coeff(2, 6) == -1);
      REQUIRE(jac_matrix->coeff(3, 7) == -1);
      REQUIRE(jac_matrix->coeff(0, 9) == -1);
      REQUIRE(jac_matrix->coeff(0, 11) == -1);
      REQUIRE(jac_matrix->coeff(0, 8) == 1);
      REQUIRE(jac_matrix->coeff(1, 10) == -1);
      REQUIRE(jac_matrix->coeff(1, 9) == 1);
      REQUIRE(jac_matrix->coeff(2, 10) == 1);
      REQUIRE(jac_matrix->coeff(2, 11) == 1);
      REQUIRE(jac_matrix->coeff(3, 8) == -1);
      REQUIRE(jac_matrix->coeff(1, 12) == -1);
      REQUIRE(jac_matrix->coeff(4, 0) == 0);
      REQUIRE(jac_matrix->coeff(5, 1) == 0);
      REQUIRE(jac_matrix->coeff(6, 2) == 0);
      REQUIRE(jac_matrix->coeff(7, 3) == 1);
      REQUIRE(jac_matrix->coeff(4, 4) == 1);
      REQUIRE(jac_matrix->coeff(5, 5) == 1);
      REQUIRE(jac_matrix->coeff(6, 6) == 1);
      REQUIRE(jac_matrix->coeff(7, 7) == 0);
      REQUIRE(jac_matrix->coeff(8, 3) == -1);
      REQUIRE(jac_matrix->coeff(9, 0) == -1);
      REQUIRE(jac_matrix->coeff(10, 1) == -1);
      REQUIRE(jac_matrix->coeff(11, 0) == -1);
      REQUIRE(jac_matrix->coeff(8, 0) == 1);
      REQUIRE(jac_matrix->coeff(9, 1) == 1);
      REQUIRE(jac_matrix->coeff(10, 2) == 1);
      REQUIRE(jac_matrix->coeff(11, 2) == 1);
      REQUIRE(jac_matrix->coeff(12, 12) == 1);

      REQUIRE(jac_matrix->coeff(8, 8) ==
              Approx(0.5932445506024201).epsilon(tolerance));
      REQUIRE(jac_matrix->coeff(9, 9) ==
              Approx(1.576675297968889).epsilon(tolerance));
      REQUIRE(jac_matrix->coeff(10, 10) ==
              Approx(2.6889982861972443).epsilon(tolerance));
      REQUIRE(jac_matrix->coeff(11, 11) ==
              Approx(13.847781018546746).epsilon(tolerance));
      REQUIRE(jac_matrix->coeff(12, 1) == Approx(-1e-11).epsilon(tolerance));
    }

    SECTION("CASE 2 (discharge between HWQ1 and HWQ2) ") {
      double init_discharge = 0.0003;
      mesh->iterate_over_links(std::bind(
          &pipenetwork::Link::update_sim_discharge, std::placeholders::_1,
          init_discharge));  // initialze discharge
      auto assembler =
          std::make_shared<pipenetwork::MatrixAssembler>(mesh, pdd_mode);
      auto jac_matrix = assembler->jac_matrix();

      // test residual (headloss part)
      assembler->assemble_residual();
      auto res_vec = assembler->residual_vector();
      REQUIRE(res_vec->coeff(8) ==
              Approx(-88.39196304144197).epsilon(tolerance));
      REQUIRE(res_vec->coeff(9) ==
              Approx(9.822533629748455e-05).epsilon(tolerance));
      REQUIRE(res_vec->coeff(10) ==
              Approx(-3.047832478024294).epsilon(tolerance));
      REQUIRE(res_vec->coeff(11) ==
              Approx(-3.047137296722242).epsilon(tolerance));

      // test jacobians (jac_g part)
      assembler->update_jacobian();
//      std::cout << (*jac_matrix) << std::endl;
      REQUIRE(jac_matrix->coeff(8, 8) ==
              Approx(0.17061370647304652).epsilon(tolerance));
      REQUIRE(jac_matrix->coeff(9, 9) ==
              Approx(0.4534427096174827).epsilon(tolerance));
      REQUIRE(jac_matrix->coeff(10, 10) ==
              Approx(0.7733403768174626).epsilon(tolerance));
      REQUIRE(jac_matrix->coeff(11, 11) ==
              Approx(3.9825418431609636).epsilon(tolerance));
      REQUIRE(jac_matrix->coeff(12, 1) == Approx(-1e-11).epsilon(tolerance));
    }

      SECTION("CASE 3 (discharge smaller than HWQ1) ") {
          double init_discharge = 0.0001;
          mesh->iterate_over_links(std::bind(
                  &pipenetwork::Link::update_sim_discharge, std::placeholders::_1,
                  init_discharge));  // initialze discharge
          auto assembler =
                  std::make_shared<pipenetwork::MatrixAssembler>(mesh, pdd_mode);
          auto jac_matrix = assembler->jac_matrix();

          // test residual (headloss part)
          assembler->assemble_residual();
          auto res_vec = assembler->residual_vector();
          REQUIRE(res_vec->coeff(8) ==
                  Approx(-88.39198847627793).epsilon(tolerance));
          REQUIRE(res_vec->coeff(9) ==
                  Approx(3.0626775928768335e-05).epsilon(tolerance));
          REQUIRE(res_vec->coeff(10) ==
                  Approx(-3.0479477664500183).epsilon(tolerance));
          REQUIRE(res_vec->coeff(11) ==
                  Approx(-3.0477310080985545).epsilon(tolerance));

          // test jacobians (jac_g part)
          assembler->update_jacobian();
//          std::cout << (*jac_matrix) << std::endl;
          REQUIRE(jac_matrix->coeff(8, 8) ==
                  Approx(0.11523722066090049).epsilon(tolerance));
          REQUIRE(jac_matrix->coeff(9, 9) ==
                  Approx(0.30626775928768335).epsilon(tolerance));
          REQUIRE(jac_matrix->coeff(10, 10) ==
                  Approx(0.5223354998349833).epsilon(tolerance));
          REQUIRE(jac_matrix->coeff(11, 11) ==
                  Approx(2.6899190144732716).epsilon(tolerance));
          REQUIRE(jac_matrix->coeff(12, 1) == Approx(-1e-11).epsilon(tolerance));
      }


      SECTION("CASE 4 (discharge is negative) ") {
          double init_discharge = -1e-3;
          mesh->iterate_over_links(std::bind(
                  &pipenetwork::Link::update_sim_discharge, std::placeholders::_1,
                  init_discharge));  // initialze discharge
          auto assembler =
                  std::make_shared<pipenetwork::MatrixAssembler>(mesh, pdd_mode);
          auto jac_matrix = assembler->jac_matrix();

          // test residual (headloss part)
          assembler->assemble_residual();
          auto res_vec = assembler->residual_vector();
          REQUIRE(res_vec->coeff(8) ==
                  Approx(-88.3923203264312).epsilon(tolerance));
          REQUIRE(res_vec->coeff(9) ==
                  Approx(-0.000851336553978881).epsilon(tolerance));
          REQUIRE(res_vec->coeff(10) ==
                  Approx(-3.0494519429191147).epsilon(tolerance));
          REQUIRE(res_vec->coeff(11) ==
                  Approx( -3.055477203573731).epsilon(tolerance));

          // test jacobians (jac_g part)
          assembler->update_jacobian();
//          std::cout << (*jac_matrix) << std::endl;
          REQUIRE(jac_matrix->coeff(8, 8) ==
                  Approx(0.5932445506024201).epsilon(tolerance));
          REQUIRE(jac_matrix->coeff(9, 9) ==
                  Approx(1.576675297968889).epsilon(tolerance));
          REQUIRE(jac_matrix->coeff(10, 10) ==
                  Approx(2.6889982861972443).epsilon(tolerance));
          REQUIRE(jac_matrix->coeff(11, 11) ==
                  Approx(13.847781018546746).epsilon(tolerance));
          REQUIRE(jac_matrix->coeff(12, 1) == Approx(-1e-11).epsilon(tolerance));
      }
  }
}
