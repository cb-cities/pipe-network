#include "eigen_gmres.h"

// Eigen GMRES Solver
bool pipenetwork::EigenGMRES::solve() {
  bool convergence = false;
  const size_t n = vec_b_->size();
  Eigen::VectorXd x_diff(n);

  // Apply restraints (only works when the known values are 0)
  (*vec_b_) = (*vec_b_).array() * vrestraints_.array();
/*
  // Apply Drichlet boundary condition (works for arbitrary known values)
  Eigen::VectorXd x_fixed(n);
  x_fixed = (*vec_x_).array() - (*vec_x_).array() * vrestraints_.array();
  (*vec_b_) -= (*mat_a_) * x_fixed;
  for (int i = 0; i < n; i++) {
    if (vrestraints_(i) == 0) {
      vec_b_->coeffRef(i) = x_fixed(i);
      for (int j = 0; j < n; j++) {
        mat_a_->coeffRef(i, j) = 0;
        mat_a_->coeffRef(j, i) = 0;
      }
      mat_a_->coeffRef(i, i) = 1;
    }
  }
*/
  // Use GMRES solver in Eigen to solve
  Eigen::GMRES<Eigen::SparseMatrix<double>> gmres(*mat_a_);
  gmres.setMaxIterations(max_iter_);
  gmres.setTolerance(tolerance_);

  x_diff = gmres.solve(*vec_b_);
  x_diff = x_diff.array() * vrestraints_.array();

  if (gmres.info() == 0) {
    convergence = true;

    // Update x (variables) in the system
    (*vec_x_) += x_diff;

    // Update delta and iterations
    this->delta_ = gmres.error();
    this->niterations_ = gmres.iterations();

    std::cout << "Iteration: " << niterations_ << " of " << max_iter_
              << "; delta: " << delta_ << std::endl;
  } else if (gmres.info() == 1) {
    throw std::runtime_error(
        "The provided data did not satisfy the prerequisites.");
  } else if (gmres.info() == 2) {
    throw std::runtime_error("Iterative procedure did not converge.");
  } else {
    throw std::runtime_error(
        "The inputs are invalid, or the algorithm has been improperly called.");
  }

  return convergence;
}
