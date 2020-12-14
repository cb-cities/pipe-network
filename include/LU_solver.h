#ifndef PIPE_NETWORK_MKL_UNSYM_H
#define PIPE_NETWORK_MKL_UNSYM_H

#include <Eigen/Sparse>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>

#include "mkl_pardiso.h"
#include "solver.h"

namespace pipenetwork {
namespace linear_system {

//! Pipe network LU Solver class
//! \brief Linear system solver using Eigen LU solver
class LUSolver : public Solver {
 public:
  //! default constructor
  LUSolver() = default;
  //! Call MKL Pardiso Unsymmetric solver
  //! \retval status Return status of the solver
  Eigen::VectorXd solve() override;

 protected:
  //! Return the type of assembler for Factory
  //! \retval assembler_type Eigen assembler
  std::string assembler_type() const override { return "LU Solver"; }

 private:
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      LU_solver_;
};
}  // namespace linear_system
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_MKL_UNSYM_H
