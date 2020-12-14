#include "LU_solver.h"
#include "factory.h"

static Register<pipenetwork::linear_system::Solver,
                pipenetwork::linear_system::LUSolver>
    registry("LU");

Eigen::VectorXd pipenetwork::linear_system::LUSolver::solve() {
  //    std::cout << "residual in linear solver: " << *vec_b_ << std::endl;
  LU_solver_.analyzePattern(matrix_assembler_->jac_matrix());
  // Compute the numerical factorization
  LU_solver_.factorize(matrix_assembler_->jac_matrix());
  // Use the factors to solve the linear system
  return LU_solver_.solve(matrix_assembler_->residual_vector());
}