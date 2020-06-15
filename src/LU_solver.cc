#include "LU_solver.h"
#include "factory.h"

static Register<pipenetwork::linear_system::Solver,
                pipenetwork::linear_system::LUSolver>
    registry("LU");

Eigen::VectorXd pipenetwork::linear_system::LUSolver::solve() {
  LU_solver_.analyzePattern(*mat_a_);
  // Compute the numerical factorization
  LU_solver_.factorize(*mat_a_);
  // Use the factors to solve the linear system
  return LU_solver_.solve(*vec_b_);
}