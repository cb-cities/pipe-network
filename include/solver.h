#ifndef PIPE_NETWORK_SOLVER_H_
#define PIPE_NETWORK_SOLVER_H_

#include "matrix_assembler.h"
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
      const std::shared_ptr<linear_system::MatrixAssembler>& matrix_assember) {
    matrix_assembler_ = matrix_assember;

    // pointers for sparse matrix info
    ia_ = matrix_assembler_->jac_matrix().outerIndexPtr();
    ja_ = matrix_assembler_->jac_matrix().innerIndexPtr();
    a_ = matrix_assembler_->jac_matrix().valuePtr();

    // matrix a basic information
    rowsA_ = matrix_assembler_->jac_matrix().rows();
    colsA_ = matrix_assembler_->jac_matrix().cols();
    nnzA_ = matrix_assembler_->jac_matrix().nonZeros();
  }

  //! Solve
  virtual Eigen::VectorXd solve() = 0;

  virtual std::string assembler_type() const = 0;

 protected:
  std::shared_ptr<linear_system::MatrixAssembler> matrix_assembler_;

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
