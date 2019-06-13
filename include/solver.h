#ifndef PIPE_NETWORK_SOLVER_H_
#define PIPE_NETWORK_SOLVER_H_

#include <Eigen/Sparse>


namespace pipenetwork {

//! Pipe network solver base class
//! \brief Base class for the pipe network solver
class Solver {
 public:
  //! Constructor with max iterations and tolerance
  Solver(unsigned max_iter, double tolerance)
      : max_iter_{max_iter}, tolerance_{tolerance} {}

  // using VectorX = typename Eigen::VectorXd;
  // using SparseMatrix = typename Eigen::SparseMatrix<double>;

  //! Copy Matrix A, Vectors b and x
  //! \param[in] mat_a_ pointer to assembled A matrix
  //! \param[in] vec_x_ pointer to assembled x vector
  //! \param[in] vec_b_ pointer to assembled b vector
  void assembled_matrices(std::shared_ptr<Eigen::SparseMatrix<double>> mat_a,
                          std::shared_ptr<Eigen::VectorXd> vec_x,
                          std::shared_ptr<Eigen::VectorXd> vec_b) {
    mat_a_ = mat_a;
    vec_x_ = vec_x;
    vec_b_ = vec_b;
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
  unsigned max_iter_;
  //! Tolerance
  double tolerance_;
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
  std::shared_ptr<Eigen::SparseMatrix<double>> mat_a_;
};
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_SOLVER_H_
