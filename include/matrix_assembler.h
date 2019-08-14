
#ifndef PIPE_NETWORK_MATRIX_ASSEMBLER_H
#define PIPE_NETWORK_MATRIX_ASSEMBLER_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>

#include "curves.h"
#include "mesh.h"
#include "settings.h"

namespace pipenetwork {

//! Mtrix assembler class
//! \brief Class for assembling matrix for solving
class MatrixAssembler {

 public:
  //! Constructor
  explicit MatrixAssembler(const std::shared_ptr<Mesh>& mesh,
                           std::shared_ptr<Curves>& curves_info,
                           bool pdd_mode = false)
      : mesh_{mesh}, curves_info_{curves_info}, pdd_{pdd_mode} {
    nnodes_ = mesh_->nnodes();
    nlinks_ = mesh_->nlinks();
    npumps_ = mesh_->npumps();
    npipes_ = mesh_->npipes();
    nvalves_ = mesh_->nvalves();

    init_variable_vector();
    assemble_balance_headloss_matrix();
    initialize_jacobian();
  }

  //! Destructor
  ~MatrixAssembler() = default;

  //! get the variable vector
  std::shared_ptr<Eigen::VectorXd> variable_vector() const {
    return variable_vec_;
  }

  //! Method to assemble residual from the variable vector
  void assemble_residual();
  std::shared_ptr<Eigen::VectorXd> residual_vector() const {
    return residual_vec_;
  }
  //! Method to update jacobian matrix from the variable vector
  void update_jacobian();
  std::shared_ptr<Eigen::SparseMatrix<double, Eigen::RowMajor>> jac_matrix()
      const {
    return jac_;
  }

 private:
  //! the mesh ptr
  std::shared_ptr<Mesh> mesh_;
  //! the curves info ptr
  std::shared_ptr<Curves> curves_info_;
  //! basic info from mesh for matrix assembling
  Index nnodes_{0};
  Index nlinks_{0};
  Index npipes_{0};
  Index npumps_{0};
  Index nvalves_{0};

  //! if it is pressure demand simulation mode
  bool pdd_{false};
  //! nodal id to corresponding matrix entry number (0-nnodes_)
  std::map<std::string, Index> node_id_map_;
  //! link id to corresponding matrix entry number (0-nlinks_)
  std::map<std::string, Index> link_id_map_;

  //! demand (for junction) and heads (for sources) of nodes
  Eigen::VectorXd demands_heads_vec_;
  //! junction elevations vector
  Eigen::VectorXd elevations_;
  //! matrix entry id for node sources (reservoir/tank)
  std::vector<Index> source_idx_;
  //! nodal ids for possible leak nodes
  std::vector<std::string> leak_ids_;
  std::vector<double> leak_area_;
  //! resistence coefficients for links
  Eigen::VectorXd link_resistance_coeff_vec_;
  //! minor loss coefficients for links
  Eigen::VectorXd link_minor_loss_coeff_vec_;

  //! node balance matrix
  Eigen::SparseMatrix<double> node_balance_mat_;
  //! link headloss matrix
  Eigen::SparseMatrix<double> headloss_mat_;
  //! internal connectivity graph
  Eigen::SparseMatrix<int> internal_graph_;

  //! map of for sub-jacobian triplets (row, col, value)
  std::map<std::string, std::vector<Eigen::Triplet<double>>> sub_jac_trip_;

  //! variable vector
  std::shared_ptr<Eigen::VectorXd> variable_vec_{
      std::make_shared<Eigen::VectorXd>()};
  //! residual vector
  std::shared_ptr<Eigen::VectorXd> residual_vec_{
      std::make_shared<Eigen::VectorXd>()};
  //! Jacobian matrix
  std::shared_ptr<Eigen::SparseMatrix<double, Eigen::RowMajor>> jac_{
      std::make_shared<Eigen::SparseMatrix<double, Eigen::RowMajor>>()};

  //! Initialize variable vector
  void init_variable_vector();

  //! Initialize matrices that contain information about node balance and link
  //! headloss (include sub_jacobian_b and sub_jacobian_f)
  void assemble_balance_headloss_matrix();

  //! Assemble residual for leakage equation
  void assemble_leak_residual();
  //! Assemble residual for mass conservation equation
  void assemble_demand_head_residual();
  //! Assemble residual for energy conservation equation (pipes)
  void assemble_headloss_residual_pipe();
  //! Assemble residual for energy conservation equation (pumps)
  void assemble_headloss_residual_pump();
  //! Assemble residual for energy conservation equation (valves)
  void assemble_headloss_residual_valve();

  //! Method to Set the jacobian entries that depend on the network status but
  //! do not depend on the value of any variable. (status of valves etc.)
  void set_jac_const();

  //! Method to update jacobian d part (pressure-demand equation)
  void update_jac_d();
  //! Method to update jacobian f part (for power pumps only)
  void update_jac_f();
  //! Method to update jacobian g part (harzen-williams headloss equation)
  void update_jac_g_pipe();
  void update_jac_g_pump();
  void update_jac_g_valve();
  //! Method to update jacobian h part (leakage equation)
  void update_jac_h();

  //! Assemble Jacobian matrix for nodal head, demand and pipe discharge as
  //! variables
  //!    Create the jacobian as a sparse matrix
  //!           Initialize all jacobian entries that have the possibility to be
  //!           non-zero
  //!            Structure of jacobian:
  //!            H_n => Head for node id n
  //!            D_n => Demand for node id n
  //!            F_l => Flow for link id l
  //!            node_bal_n => node balance for node id n
  //!            D/H_n      => demand/head equation for node id n
  //!            headloss_l => headloss equation for link id l
  //!            in link refers to a link that has node_n as an end node
  //!    out link refers to a link that has node_n as a start node
  //!            Note that there will only be leak variables and equations for
  //!            nodes with leaks. Thus some of the rows and columns below may
  //!            be missing. The leak id is equal to the node id though.
  //!    Variable          H_1   H_2   H_n   H_(N-1)   H_N   D_1   D_2   D_n
  //!    D_(N-1)   D_N   F_1   F_2   F_l   F_(L-1)   F_L      Dleak1  Dleak2
  //!    Dleakn  Dleak(N-1)  DleakN
  //!            Equation
  //!    node_bal_1         0     0     0     0         0     -1    0     0 0 0
  //!    (1 for in link, -1 for out link)       -1      0        0        0 0
  //!    node_bal_2         0     0     0     0         0     0     -1    0 0 0
  //!    (1 for in link, -1 for out link)        0     -1        0        0 0
  //!    node_bal_n         0     0     0     0         0     0     0     -1 0
  //!    0    (1 for in link, -1 for out link)        0      0       -1        0
  //!    0 node_bal_(N-1)     0     0     0     0         0     0     0     0 -1
  //!    0    (1 for in link, -1 for out link)        0      0        0       -1
  //!    0 node_bal_N         0     0     0     0         0     0     0     0 0
  //!    -1   (1 for in link, -1 for out link)        0      0        0        0
  //!    -1 D/H_1              *1    0     0     0         0     *2    0     0
  //!    0         0     0      0     0    0         0          0      0 0 0 0
  //!    D/H_2              0     *1    0     0         0     0     *2    0 0 0
  //!    0      0     0    0         0          0      0        0        0 0
  //!    D/H_n              0     0     *1    0         0     0     0     *2 0
  //!    0     0      0     0    0         0          0      0        0        0
  //!    0 D/H_(N-1)          0     0     0     *1        0     0     0     0 *2
  //!    0     0      0     0    0         0          0      0        0        0
  //!    0 D/H_N              0     0     0     0         *1    0     0     0 0
  //!    *2    0      0     0    0         0          0      0        0        0
  //!    0 headloss_1         (NZ for start/end node *3    )    0     0     0 0
  //!    0     *4     0     0    0         0          0      0        0        0
  //!    0 headloss_2         (NZ for start/end node *3    )    0     0     0 0
  //!    0     0      *4    0    0         0          0      0        0        0
  //!    0 headloss_l         (NZ for start/end node *3    )    0     0     0 0
  //!    0     0      0     *4   0         0          0      0        0        0
  //!    0 headloss_(L-1)     (NZ for start/end node *3    )    0     0     0 0
  //!    0     0      0     0    *4        0          0      0        0        0
  //!    0 headloss_L         (NZ for start/end node *3    )    0     0     0 0
  //!    0     0      0     0    0         *4         0      0        0        0
  //!    0 leak flow 1        *5    0     0     0         0     0     0     0 0
  //!    0     0      0     0    0         0          1      0        0        0
  //!    0 leak flow 2        0     *5    0     0         0     0     0     0 0
  //!    0     0      0     0    0         0          0      1        0        0
  //!    0 leak flow n        0     0     *5    0         0     0     0     0 0
  //!    0     0      0     0    0         0          0      0        1        0
  //!    0 leak flow N-1      0     0     0     *5        0     0     0     0 0
  //!    0     0      0     0    0         0          0      0        0        1
  //!    0 leak flow N        0     0     0     0         *5    0     0     0 0
  //!    0     0      0     0    0         0          0      0        0        0
  //!    1 *1: 1 for tanks and reservoirs 1 for isolated junctions 0 for
  //!    junctions if the simulation is demand-driven and the junction is not
  //!    isolated f(H) for junctions if the simulation is pressure dependent
  //!    demand and the junction is not isolated *2: 0 for tanks and reservoirs
  //!    1 for non-isolated junctions
  //!    0 for isolated junctions
  //!    *3: 0 for closed/isolated links
  //!    pipes   head_pumps  power_pumps  active_PRV   open_prv active/openTCV
  //!    active_FCV   open_FCV
  //!            start node    -1        1            f(F)        0 -1 -1 0 -1
  //!    end node       1       -1            f(F)        1              1 1 0 1
  //!    *4: 1 for closed/isolated links
  //!    f(F) for pipes
  //!    f(F) for head pumps
  //!    f(Hstart,Hend) for power pumps
  //!    0 for active PRVs
  //!    f(F) for open PRVs
  //!    f(F) for open or active TCVs
  //!    f(F) for open FCVs
  //!    1 for active FCVs
  //!    *5: 0 for inactive leaks
  //!    0 for leaks at isolated junctions
  //!    f(H-z) otherwise

  void initialize_jacobian();
};

}  // namespace pipenetwork
#endif  // PIPE_NETWORK_MATRIX_ASSEMBLER_H
