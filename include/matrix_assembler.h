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

  //! Obtain global nodal and pipe indices and pointers from meshes
  //! \param[in] mesh meshes that form the pipe network
  void global_nodal_pipe_indices(const std::shared_ptr<Mesh>& mesh);

  //! Return number of nodes in the network
  //! \retval nnode_ number of nodes in the network
  unsigned nnodes() { return nnode_; }

  //! Return number of pipes in the network
  //! \retval nnode_ number of pipes in the network
  unsigned npipes() { return npipe_; }

  //! Assign pipe roughness coefficient for Hazen-Williams equation
  //! \param[in] mesh mesh where the pipes were created
  //! \param[in] pipe_roughness vector of pair of pipe index and roughness
  void assign_pipe_roughness(
      const std::shared_ptr<Mesh>& mesh,
      const std::vector<std::pair<Index, double>>& pipe_roughness);

  //! Initialize discharges in pipes
  //! \param[in] mesh mesh where the pipes were created
  void initialize_pipe_discharge(const std::shared_ptr<Mesh>& mesh);

  //! Assign initial heads for nodes that have known head
  //! \param[in] mesh mesh where the nodes were created
  //! \param[in] node_head vector of pair of nodal index and initial nodal head
  void assign_node_head(const std::shared_ptr<Mesh>& mesh,
                        const std::vector<std::pair<Index, double>>& node_head);

  //! Assign initial discharges for nodes that have known discharge
  //! \param[in] mesh mesh where the nodes were created
  //! \param[in] node_discharge vector of pair of nodal index and initial
  //! discharge
  void assign_node_discharge(
      const std::shared_ptr<Mesh>& mesh,
      const std::vector<std::pair<Index, double>>& node_discharge);

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
  //! \retval node_head_vec_ pointer to nodal head vector
  std::shared_ptr<Eigen::VectorXd> node_head_vec() const {
    return node_head_vec_;
  }

  //! Return nodal discharge vector
  //! \retval node_discharge_vec_ pointer to nodal discharge vector
  std::shared_ptr<Eigen::VectorXd> node_discharge_vec() const {
    return node_discharge_vec_;
  }

  //! Return pipe discharge vector
  //! \retval pipe_discharge_vec_ pointer to pipe discharge vector
  std::shared_ptr<Eigen::VectorXd> pipe_discharge_vec() const {
    return pipe_discharge_vec_;
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
