#include "eigen_cg_lem.h"

// Eigen CG Solver
bool pipenetwork::EigenCG_LEM::solve() {
  bool convergence = false;
  double init_delta, new_delta;

  const size_t n = vec_b_->size();
  Eigen::VectorXd r(n), p(n), ap(n),x(n);
  x.setZero();

  r = ( (*vec_b_) - (*mat_a_) * x );
  p = r;
  auto r_k_norm = r.dot(r);

  std::cout << "init_delta " << r_k_norm << std::endl;
    std::cout << (*mat_a_) * p << std::endl;

  unsigned i = 0;
  double alpha = 0., r_kplus1_norm = 0., beta = 0.;
  for (i = 0; i < 1; ++i) {
    ap = (*mat_a_) * p;
    alpha = r_k_norm / p.dot(ap);
    std::cout<<ap<<std::endl;

    x += alpha * p;
    r -= alpha * ap;

    r_kplus1_norm = r.dot(r);
    std::cout << i << " " << r_k_norm<< " "<<r_kplus1_norm << std::endl;
    beta = r_kplus1_norm / r_k_norm;
    r_k_norm = r_kplus1_norm;

    // Break if convergence criterion is achieved
    if (r_k_norm < tolerance_) {
      convergence = true;
      break;
    }
    p = r+beta*p;
    std::cout << i << " " << r_kplus1_norm << std::endl;
  }

  // Update delta and iterations
  this->delta_ = r_k_norm;
  this->niterations_ = i;

  std::cout << i << " " << r_k_norm << std::endl;
  return convergence;
}

//! Precondition Jacobian

Eigen::VectorXd pipenetwork::EigenCG_LEM::precondition_jacobian() {
  const size_t n = vec_b_->size();
  Eigen::VectorXd vm(n);

  vm.setZero();
  for (size_t i = 0; i < n; ++i) {
    vm[i] += (*mat_a_).coeff(i, i);
  }
  for (size_t i = 0; i < n; ++i) {
    vm[i] = 1 / vm[i];
    // When beta is zero, the stiffness matrix will have zero value elements
    if (!std::isfinite(vm[i])) vm[i] = 1.0;
  }
  return vm;
}