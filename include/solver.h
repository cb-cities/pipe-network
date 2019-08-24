#ifndef PIPE_NETWORK_SOLVER_H_
#define PIPE_NETWORK_SOLVER_H_

#include <Eigen/Sparse>
#include <memory>

namespace pipenetwork {

//! Pipe network solver base class
//! \brief Base class for the pipe network solver
class Solver {
 public:
  //! Default constructor for noniterative solver
  Solver() = default;
  //! Constructor with max iterations and tolerance for iterative solvers
  Solver(unsigned max_iter, double tolerance)
      : max_iter_{max_iter}, tolerance_{tolerance} {}

  //! Copy Matrix A, Vectors b and x
  //! \param[in] mat_a_ pointer to assembled A matrix
  //! \param[in] vec_x_ pointer to assembled x vector
  //! \param[in] vec_b_ pointer to assembled b vector
  void assembled_matrices(
      std::shared_ptr<Eigen::SparseMatrix<double, Eigen::RowMajor>> mat_a,
      std::shared_ptr<Eigen::VectorXd> vec_x,
      std::shared_ptr<Eigen::VectorXd> vec_b) {
    mat_a_ = mat_a;
    vec_x_ = vec_x;
    vec_b_ = vec_b;

    ia_ = mat_a_->outerIndexPtr();
    ja_ = mat_a_->innerIndexPtr();
    a_ = mat_a_->valuePtr();
  }
  //! Restrain vector
  void restrains(const Eigen::VectorXd& restraints) {
    vrestraints_ = restraints;
  }

  //! Solve
  virtual bool solve() = 0;
  virtual std::string assembler_type() const = 0;

  //! number of iterations
  unsigned niterations() const { return niterations_; }

  //! Delta
  double delta() const { return delta_; }

 protected:
  //! Maximum number of iterations
  unsigned max_iter_{1000};
  //! Tolerance
  double tolerance_{1E-6};
  //! Displacement restraint
  Eigen::VectorXd vrestraints_;
  //! Iterations
  unsigned niterations_{0};
  //! Delta
  double delta_{std::numeric_limits<double>::quiet_NaN()};
  //! Displacement vector
  std::shared_ptr<Eigen::VectorXd> vec_x_;
  //! Force vector
  std::shared_ptr<Eigen::VectorXd> vec_b_;
  //! Sparse Stiffness Matrix
  std::shared_ptr<Eigen::SparseMatrix<double, Eigen::RowMajor>> mat_a_;

  //! Row index for csr format matrix
  int* ia_{NULL};
  //! Column index for csr format matrix
  int* ja_{NULL};
  //! Matrix values for csr format matrix
  double* a_{NULL};
};
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_SOLVER_H_
