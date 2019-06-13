#include "eigen_gmres.h"

// Eigen GMRES Solver
bool pipenetwork::EigenGMRES::solve() {
  bool convergence = false;
  const size_t n = vec_b_->size();
  Eigen::VectorXd x_diff(n);

  // Use GMRES solver in Eigen to solve
  Eigen::GMRES<Eigen::SparseMatrix<double>> gmres(*mat_a_);
  gmres.setMaxIterations(max_iter_);
  gmres.setTolerance(tolerance_);

  x_diff = gmres.solve(*vec_b_);
  x_diff = x_diff.array();

  if (gmres.info() == 0) {
    convergence = true;
//    std::cout << "difference :" <<x_diff<<std::endl;

    // Update x (variables) in the system
    (*vec_x_) -= x_diff;

    // Update delta and iterations
    this->delta_ = gmres.error();
    this->niterations_ = gmres.iterations();

//    std::cout << "Iteration: " << niterations_ << " of " << max_iter_
//              << "; delta: " << delta_ << std::endl;
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
