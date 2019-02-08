#include "catch.hpp"

#include "eigen_gmres.h"
#include "matrix_assembler.h"
#include "settings.h"
#include "solver.h"

// Check Eigen GMRES class (solve Ax=b)
TEST_CASE("Eigen GMRES solver is checked", "[EigenGMRES]") {
  // Tolerance
  const double tolerance = 1.e-10;

  // Max number of iteration
  const unsigned max_iter = 5000;

  // Input vector b
  auto b = std::make_shared<Eigen::VectorXd>();
  b->resize(8);
  b->coeffRef(0) = 3.0;
  b->coeffRef(1) = 5.1;
  b->coeffRef(2) = 22.2;
  b->coeffRef(3) = 13.0;
  b->coeffRef(4) = 0.0;
  b->coeffRef(5) = 140.4;
  b->coeffRef(6) = 50.5;
  b->coeffRef(7) = 27.0;

  // Input initial guess vector x
  auto x = std::make_shared<Eigen::VectorXd>();
  x->resize(8);
  x->setZero();

  // Input matrix A
  auto A = std::make_shared<Eigen::SparseMatrix<double>>();
  std::vector<Eigen::Triplet<double>> update;
  update.emplace_back(0, 0, 11.1);
  update.emplace_back(0, 1, 12.2);
  update.emplace_back(0, 3, 14.3);
  update.emplace_back(1, 1, 22.4);
  update.emplace_back(1, 2, 23.5);
  update.emplace_back(1, 4, 25.6);
  update.emplace_back(2, 0, 31.7);
  update.emplace_back(2, 2, 33.8);
  update.emplace_back(2, 3, 34.9);
  update.emplace_back(3, 1, 42.0);
  update.emplace_back(3, 4, 45.1);
  update.emplace_back(3, 5, 46.2);
  update.emplace_back(4, 4, 55.3);
  update.emplace_back(5, 4, 65.4);
  update.emplace_back(5, 5, 66.5);
  update.emplace_back(5, 6, 67.6);
  update.emplace_back(6, 4, 75.7);
  update.emplace_back(6, 6, 77.8);
  update.emplace_back(6, 7, 78.9);
  update.emplace_back(7, 6, 87.0);
  update.emplace_back(7, 7, 88.1);
  A->resize(8, 8);
  A->setFromTriplets(update.begin(), update.end());

  // Creat a eigen gmres solver and solve
  double gmres_tolerance = 1.e-25;
  auto solver =
      std::make_shared<pipenetwork::EigenGMRES>(max_iter, gmres_tolerance);
  solver->assembled_matrices(A, x, b);

  // Check the case that there is no known values in x
  SECTION("Check the case that there is no known values in x") {

    // Apply restraints and solve
    Eigen::VectorXd restraints(8);
    restraints << 1, 1, 1, 1, 1, 1, 1, 1;
    solver->restrains(restraints);
    bool issolved = solver->solve();

    // Check results
    REQUIRE(x->coeff(0) == Approx(-3470.91931206031813).epsilon(tolerance));
    REQUIRE(x->coeff(1) == Approx(-258.21997602699554).epsilon(tolerance));
    REQUIRE(x->coeff(2) == Approx(246.35010480871063).epsilon(tolerance));
    REQUIRE(x->coeff(3) == Approx(2914.71944555236860).epsilon(tolerance));
    REQUIRE(x->coeff(4) == Approx(0.0).epsilon(tolerance));
    REQUIRE(x->coeff(5) == Approx(235.02681803319942).epsilon(tolerance));
    REQUIRE(x->coeff(6) == Approx(-229.12549407112076).epsilon(tolerance));
    REQUIRE(x->coeff(7) == Approx(226.57114624503413).epsilon(tolerance));
  }

  // Check the case that there are known zeros in x
  SECTION("Check the case that there are known zeros in x") {

    // Apply restraints and solve
    Eigen::VectorXd restraints(8);
    restraints << 1, 1, 1, 1, 0, 1, 1, 1;
    solver->restrains(restraints);
    bool issolved = solver->solve();

    // Check results
    REQUIRE(x->coeff(0) == Approx(-3470.91931206031813).epsilon(tolerance));
    REQUIRE(x->coeff(1) == Approx(-258.21997602699554).epsilon(tolerance));
    REQUIRE(x->coeff(2) == Approx(246.35010480871063).epsilon(tolerance));
    REQUIRE(x->coeff(3) == Approx(2914.71944555236860).epsilon(tolerance));
    REQUIRE(x->coeff(4) == Approx(0.0).epsilon(tolerance));
    REQUIRE(x->coeff(5) == Approx(235.02681803319942).epsilon(tolerance));
    REQUIRE(x->coeff(6) == Approx(-229.12549407112076).epsilon(tolerance));
    REQUIRE(x->coeff(7) == Approx(226.57114624503413).epsilon(tolerance));
  }
}
