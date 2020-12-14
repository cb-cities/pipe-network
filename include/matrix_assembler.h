#ifndef PIPE_NETWORK_MATRIX_ASSEMBLER_H
#define PIPE_NETWORK_MATRIX_ASSEMBLER_H

#include "matrix_assembler_components.h"

namespace pipenetwork {
namespace linear_system {
//! Mtrix assembler class
//! \brief Class for assembling matrix for solving
class MatrixAssembler {

 public:
  //! Constructor
  //! \param[in] mesh the mesh pointer
  //! \param[in] curves_info the curve information pointer
  //! \param[in] pdd_mode if simulation type is pressure demand driven or
  //! demand driven
  MatrixAssembler(const std::shared_ptr<Mesh>& mesh,
                  std::shared_ptr<Curves>& curves_info, bool pdd_mode = false);

  //! Destructor
  ~MatrixAssembler() = default;

  //! get the variable vector
  Eigen::VectorXd& variable_vector() { return vars_->variables_vec(); }
  //! Method to assemble residual from the variable vector
  void assemble_residual();
  //! Method to get residual vector
  const Eigen::VectorXd& residual_vector() const {
    return res_->residual_vec();
  }
  //! Method to update jacobian matrix from the variable vector
  void update_jacobian();
  //! Method to get jacobian matrix
  const Eigen::SparseMatrix<double, Eigen::RowMajor>& jac_matrix() const {
    return jac_->jac_matrix();
  }

  //! Method for system update after variable changes
  void system_update() {
    assemble_residual();
    update_jacobian();
  }

 private:
  //! the mesh ptr
  std::shared_ptr<Mesh> mesh_;
  //! the curves info ptr
  std::shared_ptr<Curves> curves_info_;
  //! if it is pressure demand simulation mode
  bool pdd_{false};
  //! Variable object
  std::shared_ptr<Variables> vars_;
  //! Residual object
  std::shared_ptr<Residuals> res_;
  //! Jacobian object
  std::shared_ptr<Jacobian> jac_;
  //! Update curves info
  void update_curves_info_();
};
}  // namespace linear_system
}  // namespace pipenetwork
#endif  // PIPE_NETWORK_MATRIX_ASSEMBLER_H
