#ifndef PIPE_NETWORK_MATRIX_ASSEMBLER_H_
#define PIPE_NETWORK_MATRIX_ASSEMBLER_H_

#include <cmath>

#include <array>
#include <exception>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include <Eigen/Sparse>

#include "mesh.h"
#include "node.h"
#include "pipe.h"
#include "settings.h"

namespace pipenetwork {

//! Mtrix assembler class
//! \brief Class for assembling matrix for solving
class MatrixAssembler {

 public:
  //! Constructor
  MatrixAssembler();

  //! Destructor
  ~MatrixAssembler() = default;

  //! Obtain global nodal and pipe indices and pointers from meshes
  //! \param[in] mesh meshes that form the pipe network
  void global_nodal_pipe_indices(const std::shared_ptr<Mesh>& mesh);

  //! Return number of nodes in the network
  //! \retval nnode_ number of nodes in the network
  unsigned nnodes() { return nnode_; }

  //! Return number of pipes in the network
  //! \retval nnode_ number of pipes in the network
  unsigned npipes() { return npipe_; }

  //! Initialize variable vector
  void assemble_variable_vector();

  // Apply variables (head and discharge) back to nodes and pipes
  void apply_variables();

  //! Calculate and assemble residual vector
  void assemble_residual_vector();

  //! Assemble Jacobian matrix
  void assemble_jacobian();

  //! Return variable vector
  //! \retval variable_vec_ pointer to variable vector
  std::shared_ptr<Eigen::VectorXd> variable_vec() const {
    return variable_vec_;
  }

  //! Return residual vector
  //! \retval residual_vec_ pointer to residual vector
  std::shared_ptr<Eigen::VectorXd> residual_vec() const {
    return residual_vec_;
  }

  //! Return Jacobian matrix
  //! \retval jac_ pointer to Jacobian matrix
  std::shared_ptr<Eigen::SparseMatrix<double>> jac() const { return jac_; }

  //**************************************************************
  //***** Make all nodal discharge known value, test purpose *****
  //**************************************************************

  //! Assemble Jacobian matrix
  void sim_assemble_jacobian();

  //! Apply variables (head and discharge) to nodes and pipes
  void sim_apply_variables();

  //! Initialize variable vector
  void sim_assemble_variable_vector();

 private:
  //! global nodal id and corresponding nodal pointer
  std::map<Index, std::shared_ptr<pipenetwork::Node>> global_nodes_;
  //! global pipe id and corresponding pipe pointer
  std::map<Index, std::shared_ptr<pipenetwork::Pipe>> global_pipes_;
  //! number of nodes in the network
  unsigned nnode_{0};
  //! number of pipes in the network
  unsigned npipe_{0};
  //! variable vector
  std::shared_ptr<Eigen::VectorXd> variable_vec_;
  //! residual vector
  std::shared_ptr<Eigen::VectorXd> residual_vec_;
  //! Jacobian matrix
  std::shared_ptr<Eigen::SparseMatrix<double>> jac_;
};
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_MATRIX_ASSEMBLER_H_
