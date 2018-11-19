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

//! Mtrix assembler class
//! \brief Class for assembling matrix for solving
class MatrixAssembler {

 public:
  //! Constructor
  MatrixAssembler();

  //! Destructor
  ~MatrixAssembler() = default;

  //! Copy constructor
  MatrixAssembler(const MatrixAssembler&) = default;

  //! Assignment operator
  MatrixAssembler& operator=(const MatrixAssembler&) = default;

  //! Move constructor
  MatrixAssembler(MatrixAssembler&&) = default;

  //! Assign global index to nodes in the whole network
  //! \param[in] nodes nodal pointers created in mesh class
  void assign_global_nodal_indices(
      std::vector<std::shared_ptr<pipenetwork::Node>> nodes);

  //! Assign global index to pipes in the whole network
  //! \param[in] pipes pipe pointers created in mesh class
  void assign_global_pipe_indices(
      std::vector<std::shared_ptr<pipenetwork::Pipe>> pipes);

  //! Return number of nodes in the network
  //! \retval nnode_ number of nodes in the network
  unsigned nnodes() { return nnode_; }

  //! Return number of pipes in the network
  //! \retval nnode_ number of pipes in the network
  unsigned npipes() { return npipe_; }

  //! Initialize nodal head vector
  //! If head of the ndoe is unknown (hasn't been assigned), initialize to zero
  void assemble_node_head_vector();

  //! Initialize nodal discharge vector
  //! If discharge of the ndoe is unknown (hasn't been assigned), initialize to
  //! zero
  void assemble_node_discharge_vector();

  //! Apply head to nodes
  void apply_node_head();

  //! Apply discharge to nodes
  void apply_node_discharge();

  //! Initialize pipe discharge vector
  //! Calculated according to nodal heads at two end and Hazen-Williams equation
  //! If any of the nodal head is unknown, initialize the discharge to 0.001
  void assemble_pipe_discharge_vector();

  //! Assemble Jacobian matrix
  void assemble_jacobian();

  //! Return nodal head vector
  //! \retval *node_head_vec_ nodal head vector
  Eigen::VectorXd node_head_vec() const { return *node_head_vec_; }

  //! Return nodal discharge vector
  //! \retval *node_discharge_vec_ nodal discharge vector
  Eigen::VectorXd node_discharge_vec() const { return *node_discharge_vec_; }

  //! Return pipe discharge vector
  //! \retval *pipe_discharge_vec_ pipe discharge vector
  Eigen::VectorXd pipe_discharge_vec() const { return *pipe_discharge_vec_; }

  //! Return Jacobian matrix
  //! \retval *jac_ Jacobian matrix
  Eigen::SparseMatrix<double> jac() const { return *jac_; }

 private:
  //! global nodal id and corresponding nodal pointer
  std::map<std::shared_ptr<pipenetwork::Node>, Index> global_nodes_;
  //! global pipe id and corresponding pipe pointer
  std::map<std::shared_ptr<pipenetwork::Pipe>, Index> global_pipes_;
  //! number of nodes in the network
  unsigned nnode_{0};
  //! number of pipes in the network
  unsigned npipe_{0};
  //! nodal head vector
  std::shared_ptr<Eigen::VectorXd> node_head_vec_;
  //! nodal discharge vector
  std::shared_ptr<Eigen::VectorXd> node_discharge_vec_;
  //! pipe discharge vector
  std::shared_ptr<Eigen::VectorXd> pipe_discharge_vec_;
  //! Jacobian matrix
  std::shared_ptr<Eigen::SparseMatrix<double>> jac_;
};

#endif  // PIPE_NETWORK_MATRIX_ASSEMBLER_H_
