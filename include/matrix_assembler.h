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
  explicit MatrixAssembler(bool pdd_mode = false) ;

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

  //! Assemble Jacobian matrix for nodal head, demand and pipe discharge as
  //! variables
  void assemble_jacobian();

  //! Apply variables (head, demand and discharge) to nodes and pipes
  void apply_variables();

  //! Initialize variable vector (unknown variables only)
  void assemble_residual_vector();

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

 private:
  //! global nodal id and corresponding nodal pointer
  std::map<Index, std::shared_ptr<pipenetwork::Node>> global_nodes_;
  //! global pipe id and corresponding pipe pointer
  std::map<Index, std::shared_ptr<pipenetwork::Pipe>> global_pipes_;
  //! number of nodes in the network
  unsigned nnode_{0};
  //! number of known nodes in the network
  unsigned nnode_known_{0};
  //! number of pipes in the network
  unsigned npipe_{0};
  //! pdd mode
  bool pdd_{false};
  //! variable vector
  std::shared_ptr<Eigen::VectorXd> variable_vec_;
  //! residual vector
  std::shared_ptr<Eigen::VectorXd> residual_vec_;
  //! Jacobian matrix
  std::shared_ptr<Eigen::SparseMatrix<double>> jac_;
  //! id map that maps node id to variable id
  std::shared_ptr<std::map<Index, Index>> id_map_;

  //! Assemble pressure demand part of the jacobian matrix
  //! \param[in] n the corresponding node pointer
  //! \param[in] index the position of the corresponding node in variable vector
  //! \param[in] update the vector that stores the position of element in
  //! jacobian matrix
  void construct_demand_jac(const std::shared_ptr<pipenetwork::Node> & n, Index index,
                            std::vector<Eigen::Triplet<double>>& update);

  //! Method to get the corresponding value of jacobian for head-pressure
  //! equations \param[in] node pointer for the desired node \retval The value
  //! for the corresponding jacobian entry
  double get_pressure_head_jacob(const std::shared_ptr<pipenetwork::Node> & node);

  //! Method to assemble residuals for demand equation in pressure-demand mode
  //! \param[in] node pointer for the desired node
  //! \param[in]  ndex the position of the corresponding node in variable vector
  void assemble_pdd_residual(const std::shared_ptr<pipenetwork::Node> & node,
                             Index index);
};
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_MATRIX_ASSEMBLER_H_
