#ifndef PIPE_NETWORK_SOLVER_H_
#define PIPE_NETWORK_SOLVER_H_

#include <Eigen/Sparse>
#include <memory>

namespace pipenetwork {
namespace linear_system {

//! Pipe network solver base class
//! \brief Base class for the pipe network solver
class Solver {
 public:
  //! Default constructor for noniterative solver
  Solver() = default;

  //! Copy Matrix A, Vectors b and x
  //! \param[in] mat_a_ pointer to assembled A matrix
  //! \param[in] vec_x_ pointer to assembled x vector
  //! \param[in] vec_b_ pointer to assembled b vector
  void assembled_matrices(
      const Eigen::SparseMatrix<double, Eigen::RowMajor>& mat_a,
      Eigen::VectorXd& vec_x, const Eigen::VectorXd& vec_b) {
    mat_a_ =
        std::make_shared<Eigen::SparseMatrix<double, Eigen::RowMajor>>(mat_a);
    vec_x_ = std::make_shared<Eigen::VectorXd>(vec_x);
    vec_b_ = std::make_shared<Eigen::VectorXd>(vec_b);
    // pointers for sparse matrix info
    ia_ = mat_a_->outerIndexPtr();
    ja_ = mat_a_->innerIndexPtr();
    a_ = mat_a_->valuePtr();
    // matrix a basic information
    rowsA_ = mat_a_->rows();
    colsA_ = mat_a_->cols();
    nnzA_ = mat_a_->nonZeros();
  }

  //! Solve
  virtual Eigen::VectorXd solve() = 0;

  virtual std::string assembler_type() const = 0;

 protected:
  //! Tolerance
  double tolerance_{1E-6};
  //! Variable vector
  std::shared_ptr<Eigen::VectorXd> vec_x_;
  //! Residual vector
  std::shared_ptr<Eigen::VectorXd> vec_b_;
  //! The A Matrix
  std::shared_ptr<Eigen::SparseMatrix<double, Eigen::RowMajor>> mat_a_;

  //! Row index for csr format matrix
  int* ia_{nullptr};
  //! Column index for csr format matrix
  int* ja_{nullptr};

  int rowsA_ = 0; /* number of rows of A */
  int colsA_ = 0; /* number of columns of A */
  int nnzA_ = 0;  /* number of nonzeros of A */
  int baseA_ = 0;

  //! Matrix values for csr format matrix
  double* a_{nullptr};
};
}  // namespace linear_system
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_SOLVER_H_
