#ifndef PIPE_NETWORK_CUDA_UNSYM_H
#define PIPE_NETWORK_CUDA_UNSYM_H

#include <Eigen/Sparse>
#include <assert.h>
#include <ctype.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cusolverSp.h"
#include "cusparse.h"
#include "solver.h"

#include "helper_cuda.h"
#include "helper_cusolver.h"

namespace pipenetwork {

class Cuda_unsym : public Solver {
 public:
  //! Constructor with tolerance, precision and iterations
  //! \param[in] max_iter Maximum number of iterations
  //! \param[in] tolerance Tolerance for solver to achieve convergence
  Cuda_unsym() = default;
  //! Call Pardiso Unsymmetric solver
  //! \retval status Return status of the solver
  Eigen::VectorXd solve() override;

 protected:
  //! Vector x
  using Solver::vec_x_;
  //! Vector b
  using Solver::vec_b_;
  //! Matrix A
  using Solver::mat_a_;
  //! Row index for csr format matrix
  using Solver::ia_;
  //! Column index for csr format matrix
  using Solver::ja_;
  //! Matrix values for csr format matrix
  using Solver::a_;

  using Solver::colsA_;
  using Solver::nnzA_;
  using Solver::rowsA_;
  //! Return the type of assembler for Factory
  //! \retval assembler_type Eigen assembler
  std::string assembler_type() const override { return "Cuda"; }

 private:
  Eigen::VectorXd x_diff_;
  cusolverSpHandle_t handle_ = NULL;
  cusparseHandle_t cusparseHandle_ = NULL; /* used in residual evaluation */
  cudaStream_t stream_ = NULL;
  cusparseMatDescr_t descrA_ = NULL;

  int baseA_ = 0; /* base index in CSR format */

  /* CSR(A) from I/O */
  int* h_csrRowPtrA_ = NULL;
  int* h_csrColIndA_ = NULL;
  double* h_csrValA_ = NULL;

  double* h_z_ = NULL;  /* z = B \ (Q*b) */
  double* h_x_ = NULL;  /* x = A \ b */
  double* h_b_ = NULL;  /* b = ones(n,1) */
  double* h_Qb_ = NULL; /* Q*b */
  double* h_r_ = NULL;  /* r = b - A*x */

  int* h_Q_ = NULL; /* <int> n */
  /* reorder to reduce zero fill-in */
  /* Q = symrcm(A) or Q = symamd(A) */
  /* B = Q*A*Q' or B = A(Q,Q) by MATLAB notation */
  int* h_csrRowPtrB_ = NULL; /* <int> n+1 */
  int* h_csrColIndB_ = NULL; /* <int> nnzA */
  double* h_csrValB_ = NULL; /* <double> nnzA */
  int* h_mapBfromA_ = NULL;  /* <int> nnzA */

  size_t size_perm_ = 0;
  void* buffer_cpu_ = NULL; /* working space for permutation: B = Q*A*Q^T */

  /* device copy of A: used in residual evaluation */
  int* d_csrRowPtrA_ = NULL;
  int* d_csrColIndA_ = NULL;
  double* d_csrValA_ = NULL;

  /* device copy of B: used in B*z = Q*b */
  int* d_csrRowPtrB_ = NULL;
  int* d_csrColIndB_ = NULL;
  double* d_csrValB_ = NULL;

  int* d_Q_ = NULL;     /* device copy of h_Q */
  double* d_z_ = NULL;  /* z = B \ Q*b */
  double* d_x_ = NULL;  /* x = A \ b */
  double* d_b_ = NULL;  /* a copy of h_b */
  double* d_Qb_ = NULL; /* a copy of h_Qb */
  double* d_r_ = NULL;  /* r = b - A*x */

  double tol_ = 1.e-12;
  const int reorder_ = 1; /* use reordering */
  int singularity_ = 0;   /* -1 if A is invertible under tol. */

  /* the constants are used in residual evaluation, r = b - A*x */
  const double minus_one_ = -1.0;
  const double one_ = 1.0;

  double b_inf_ = 0.0;
  double x_inf_ = 0.0;
  double r_inf_ = 0.0;
  double A_inf_ = 0.0;
  int errors_ = 0;
  int issym_ = 0;

  void assign_mtx_ptrs();
  void init_cuda_handlers();
  void allocate_host_memories();
  void allocate_device_memories();

  void reorder_mtx_A();
  void create_mtx_B();
  void create_mtx_QB();

  void cp_hdata_ddata();
  void solve_z();
  void solve_x();
  void x_to_host();
  void release_resources();
  void release_host_resources();
  void release_device_resources();

  void create_x_diff();
};

};  // namespace pipenetwork

#endif  // PIPE_NETWORK_CUDA_UNSYM_H
