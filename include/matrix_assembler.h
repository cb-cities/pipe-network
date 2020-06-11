
#ifndef PIPE_NETWORK_MATRIX_ASSEMBLER_H
#define PIPE_NETWORK_MATRIX_ASSEMBLER_H

//#include <algorithm>

#include "curves.h"
#include "mesh.h"
#include "settings.h"

namespace pipenetwork {
namespace linear_system {

////! Jacobian class
////! \brief Class for linear system jacobian
// class Jacobian {
// public:
// private:
//};

//! Vairables class
//! \brief Class for linear system variables
class Variables {
 public:
  explicit Variables(const std::shared_ptr<Mesh>& mesh) : mesh_{mesh} {
    init_variable_vectors();
  }
  //! get the variable vector
  Eigen::VectorXd& variables_vec() { return variable_vec_; }

  //! get the variable vector
  const Eigen::VectorXd& demands_heads_vec() const {
    return demands_heads_vec_;
  }

  //! get the elevation vector
  const Eigen::VectorXd& elevations() const { return elevations_; }

  //! get the resistance loss coefficient vector
  const Eigen::VectorXd& link_resistance_coeff_vec() const {
    return link_resistance_coeff_vec_;
  }

  //! get the minor loss coefficient vector
  const Eigen::VectorXd& link_minor_loss_coeff_vec() const {
    return link_minor_loss_coeff_vec_;
  }
  //! get the leak areas for broken nodes
  const Eigen::VectorXd& leak_areas() const { return leak_areas_; }

 private:
  //! the mesh ptr
  std::shared_ptr<Mesh> mesh_;
  //! variable vector
  Eigen::VectorXd variable_vec_;
  //! demand (for junction) and heads (for sources) of nodes
  Eigen::VectorXd demands_heads_vec_;
  //! junction elevations vector
  Eigen::VectorXd elevations_;
  //! resistence coefficients for links
  Eigen::VectorXd link_resistance_coeff_vec_;
  //! minor loss coefficients for links
  Eigen::VectorXd link_minor_loss_coeff_vec_;
  //! leak areas
  Eigen::VectorXd leak_areas_;

  //! Initialize variable vector
  void init_variable_vectors();
  //! resize variables
  void resize_init_variables(Index nnodes, Index nlinks, Index nleaks);
  //! initialize node related vectors (demand_head, elevations)
  void init_nodes_vecs(const std::shared_ptr<MeshNodes>& nodes);
  //! initialize link related vectors (link resistance coefficients, minor loss
  //! coefficients)
  void init_link_vecs(const std::shared_ptr<MeshLinks>& links,
                      Index link_start_idx);
  //! initialize leak nodes information
  void init_leak_nodes(const std::shared_ptr<MeshNodes>& nodes,
                       const std::vector<Index>& leak_ids,
                       Index leak_start_idx);

  //! helper function to create minor loss and flow energy loss
  //! (harzan-williams) for links \param[in] linkmap map of link ptrs
  template <typename LinkMap>
  void update_minor_resis_coeff(const LinkMap& linkmap);

  //! helper function to calcuate flow energy loss for links
  //! \param[in] link link_ptr
  double get_link_res_coeff(const std::shared_ptr<Link>& link) { return 0; }

  //! helper function to calcuate flow energy loss for pipes
  //! \param[in] pipe pipe_ptr
  double get_link_res_coeff(const std::shared_ptr<Pipe>& pipe);

  //! helper function to calcuate flow energy loss for valves
  //! \param[in] valve valve_ptr
  double get_link_res_coeff(const std::shared_ptr<Valve>& valve);

  //! helper function to minor loss for links
  //! \param[in] link link_ptr
  double get_link_minor_coeff(const std::shared_ptr<Link>& link) { return 0; }
};

//! Residuals class
//! \brief Class for linear system residuals
class Residuals {
 public:
  Residuals(const std::shared_ptr<Mesh>& mesh,
            const std::shared_ptr<Variables>& vars,
            const std::shared_ptr<Curves>& curves);

  //! assemble the residual vector
  void assemble_residual();

  //! get the residual vector
  const Eigen::VectorXd& residual_vec() const { return residual_vec_; }

 private:
  //! the mesh ptr
  std::shared_ptr<Mesh> mesh_;
  //! curves ptr
  std::shared_ptr<Curves> curves_info_;
  //! the vriable ptr
  std::shared_ptr<Variables> vars_;
  //! node balance matrix
  Eigen::SparseMatrix<double> node_balance_mat_;
  //! link headloss matrix
  Eigen::SparseMatrix<double> headloss_mat_;
  //! Residual vector
  Eigen::VectorXd residual_vec_;
  //! Variable vector
  Eigen::VectorXd& variable_vec_;
  //! number of nodes
  Index nnodes_;
  //! number of links
  Index nlinks_;

  //! isolated junctions vector (boolean mask purpose)
  Eigen::VectorXd iso_junctions_mask_;
  //! connected junctions vector (boolean mask purpose)
  Eigen::VectorXd connect_junctions_mask_;
  //! isolated/closed links vector (boolean mask purpose)
  Eigen::VectorXd iso_links_mask_;
  //! connected links vector (boolean mask purpose)
  Eigen::VectorXd connect_links_mask_;

  //! assemble the masks for isolated nodes and links
  void assemble_iso_masks();

  //! assemble the balance headloss matrix for fast residual computation
  void assemble_balance_headloss_matrix();

  //! assemble node balance residual (0-nnodes)
  void assemble_node_balance_residual();

  //! assemble demand head residual (nnodes-2*nnodes)
  void assemble_demand_head_residual();

  //! assemble headloss residual for pipes (hazan-williams equation)
  void assemble_headloss_residual_pipe();

  //! assemble headloss residual for pumps
  void assemble_headloss_residual_pump();

  //! get pump head gain for head pumps
  double get_pump_headgain(const std::shared_ptr<Pump>& pump, double link_flow);
  //! assemble headloss residual for valves
  void assemble_headloss_residual_valve();

  //! get residual for open valves
  double get_open_valve_residual(const std::shared_ptr<Valve>& valve,
                                 double link_flow);
  //! get residual for active valves based on different valve types
  double get_active_valve_residual(const std::shared_ptr<Valve>& valve,
                                   double link_flow);
  //! assemble leak residual
  void assemble_leak_residual();
};

//
////! Mtrix assembler class
////! \brief Class for assembling matrix for solving
// class MatrixAssembler {
//
// public:
//  //! Constructor
//  //! \param[in] mesh the mesh pointer
//  //! \param[in] curves_info the curve information pointer
//  //! \param[in] pdd_mode if simulation type is pressure demand driven or
//  demand
//  //! driven
//  MatrixAssembler(const std::shared_ptr<Mesh>& mesh,
//                  std::shared_ptr<Curves>& curves_info, bool pdd_mode =
//                  false);
//
//  //! Destructor
//  ~MatrixAssembler() = default;
//
//  //! get the variable vector
//  std::shared_ptr<Eigen::VectorXd> variable_vector() { return variable_vec_; }
//
//  //! Method to assemble residual from the variable vector
//  void assemble_residual();
//  //! Method to get residual vector
//  std::shared_ptr<Eigen::VectorXd> residual_vector() const {
//    return residual_vec_;
//  }
//  //! Method to update jacobian matrix from the variable vector
//  void update_jacobian();
//  //! Method to get jacobian matrix
//  std::shared_ptr<Eigen::SparseMatrix<double, Eigen::RowMajor>> jac_matrix()
//      const {
//    return jac_;
//  }
//
// private:
//  //! the mesh ptr
//  std::shared_ptr<Mesh> mesh_;
//  //! the curves info ptr
//  std::shared_ptr<Curves> curves_info_;
//
//  //! if it is pressure demand simulation mode
//  bool pdd_{false};
//  //! demand (for junction) and heads (for sources) of nodes
//  Eigen::VectorXd demands_heads_vec_;
//  //! junction elevations vector
//  Eigen::VectorXd elevations_;
//  //! resistence coefficients for links
//  Eigen::VectorXd link_resistance_coeff_vec_;
//  //! minor loss coefficients for links
//  Eigen::VectorXd link_minor_loss_coeff_vec_;
//
//  //! node balance matrix
//  Eigen::SparseMatrix<double> node_balance_mat_;
//  //! link headloss matrix
//  Eigen::SparseMatrix<double> headloss_mat_;
//
//  //! map of for sub-jacobian triplets (row, col, value)
//  std::map<std::string, std::vector<Eigen::Triplet<double>>> sub_jac_trip_;
//
//  //! variable vector
//  std::shared_ptr<Eigen::VectorXd> variable_vec_{
//      std::make_shared<Eigen::VectorXd>()};
//  //! residual vector
//  std::shared_ptr<Eigen::VectorXd> residual_vec_{
//      std::make_shared<Eigen::VectorXd>()};
//  //! Jacobian matrix
//  std::shared_ptr<Eigen::SparseMatrix<double, Eigen::RowMajor>> jac_{
//      std::make_shared<Eigen::SparseMatrix<double, Eigen::RowMajor>>()};
//
//  //! Initialize variable vector
//  void init_variable_vectors();
//
//  //! resize variables
//  void resize_init_variables(Index nnodes, Index nlinks, Index nleaks);
//
//  //! initialize node related vectors (demand_head, elevations)
//  void init_nodes_vecs(const std::shared_ptr<MeshNodes>& nodes);
//
//  //! initialize link related vectors (
//  void init_link_vecs(const std::shared_ptr<MeshLinks>& links);
//
//  //! Initialize matrices that contain information about node balance and link
//  //! headloss (include sub_jacobian_b and sub_jacobian_f)
//  void assemble_balance_headloss_matrix();
//
//  //! Assemble residual for leakage equation
//  void assemble_leak_residual();
//  //! Assemble residual for mass conservation equation
//  void assemble_demand_head_residual();
//  //! Assemble residual for energy conservation equation (pipes)
//  void assemble_headloss_residual_pipe();
//  //! Assemble residual for energy conservation equation (pumps)
//  void assemble_headloss_residual_pump();
//  //! Assemble residual for energy conservation equation (valves)
//  void assemble_headloss_residual_valve();
//
//  //! Method to Set the jacobian entries that depend on the network status but
//  //! do not depend on the value of any variable. (status of valves etc.)
//  void set_jac_const();
//
//  //! Method to update jacobian d part (pressure-demand equation)
//  void update_jac_d();
//  //! Method to update jacobian f part (for power pumps only)
//  void update_jac_f();
//  //! Method to update jacobian g part (harzen-williams headloss equation)
//  void update_jac_g_pipe();
//  void update_jac_g_pump();
//  void update_jac_g_valve();
//  //! Method to update jacobian h part (leakage equation)
//  void update_jac_h();
//
//  //! Assemble Jacobian matrix for nodal head, demand and pipe discharge as
//  //! variables
//  //!    Create the jacobian as a sparse matrix
//  //!           Initialize all jacobian entries that have the possibility to
//  be
//  //!           non-zero
//  //!            Structure of jacobian:
//  //!            H_n => Head for node id n
//  //!            D_n => Demand for node id n
//  //!            F_l => Flow for link id l
//  //!            node_bal_n => node balance for node id n
//  //!            D/H_n      => demand/head equation for node id n
//  //!            headloss_l => headloss equation for link id l
//  //!            in link refers to a link that has node_n as an end node
//  //!    out link refers to a link that has node_n as a start node
//  //!            Note that there will only be leak variables and equations for
//  //!            nodes with leaks. Thus some of the rows and columns below may
//  //!            be missing. The leak id is equal to the node id though.
//  //!    Variable          H_1   H_2   H_n   H_(N-1)   H_N   D_1   D_2   D_n
//  //!    D_(N-1)   D_N   F_1   F_2   F_l   F_(L-1)   F_L      Dleak1  Dleak2
//  //!    Dleakn  Dleak(N-1)  DleakN
//  //!            Equation
//  //!    node_bal_1         0     0     0     0         0     -1    0     0 0
//  0
//  //!    (1 for in link, -1 for out link)       -1      0        0        0 0
//  //!    node_bal_2         0     0     0     0         0     0     -1    0 0
//  0
//  //!    (1 for in link, -1 for out link)        0     -1        0        0 0
//  //!    node_bal_n         0     0     0     0         0     0     0     -1 0
//  //!    0    (1 for in link, -1 for out link)        0      0       -1 0
//  //!    0 node_bal_(N-1)     0     0     0     0         0     0     0     0
//  -1
//  //!    0    (1 for in link, -1 for out link)        0      0        0 -1
//  //!    0 node_bal_N         0     0     0     0         0     0     0     0
//  0
//  //!    -1   (1 for in link, -1 for out link)        0      0        0 0
//  //!    -1 D/H_1              *1    0     0     0         0     *2    0     0
//  //!    0         0     0      0     0    0         0          0      0 0 0 0
//  //!    D/H_2              0     *1    0     0         0     0     *2    0 0
//  0
//  //!    0      0     0    0         0          0      0        0        0 0
//  //!    D/H_n              0     0     *1    0         0     0     0     *2 0
//  //!    0     0      0     0    0         0          0      0        0 0
//  //!    0 D/H_(N-1)          0     0     0     *1        0     0     0     0
//  *2
//  //!    0     0      0     0    0         0          0      0        0 0
//  //!    0 D/H_N              0     0     0     0         *1    0     0     0
//  0
//  //!    *2    0      0     0    0         0          0      0        0 0
//  //!    0 headloss_1         (NZ for start/end node *3    )    0     0     0
//  0
//  //!    0     *4     0     0    0         0          0      0        0 0
//  //!    0 headloss_2         (NZ for start/end node *3    )    0     0     0
//  0
//  //!    0     0      *4    0    0         0          0      0        0 0
//  //!    0 headloss_l         (NZ for start/end node *3    )    0     0     0
//  0
//  //!    0     0      0     *4   0         0          0      0        0 0
//  //!    0 headloss_(L-1)     (NZ for start/end node *3    )    0     0     0
//  0
//  //!    0     0      0     0    *4        0          0      0        0 0
//  //!    0 headloss_L         (NZ for start/end node *3    )    0     0     0
//  0
//  //!    0     0      0     0    0         *4         0      0        0 0
//  //!    0 leak flow 1        *5    0     0     0         0     0     0     0
//  0
//  //!    0     0      0     0    0         0          1      0        0 0
//  //!    0 leak flow 2        0     *5    0     0         0     0     0     0
//  0
//  //!    0     0      0     0    0         0          0      1        0 0
//  //!    0 leak flow n        0     0     *5    0         0     0     0     0
//  0
//  //!    0     0      0     0    0         0          0      0        1 0
//  //!    0 leak flow N-1      0     0     0     *5        0     0     0     0
//  0
//  //!    0     0      0     0    0         0          0      0        0 1
//  //!    0 leak flow N        0     0     0     0         *5    0     0     0
//  0
//  //!    0     0      0     0    0         0          0      0        0 0
//  //!    1 *1: 1 for tanks and reservoirs 1 for isolated junctions 0 for
//  //!    junctions if the simulation is demand-driven and the junction is not
//  //!    isolated f(H) for junctions if the simulation is pressure dependent
//  //!    demand and the junction is not isolated *2: 0 for tanks and
//  reservoirs
//  //!    1 for non-isolated junctions
//  //!    0 for isolated junctions
//  //!    *3: 0 for closed/isolated links
//  //!    pipes   head_pumps  power_pumps  active_PRV   open_prv active/openTCV
//  //!    active_FCV   open_FCV
//  //!            start node    -1        1            f(F)        0 -1 -1 0 -1
//  //!    end node       1       -1            f(F)        1              1 1 0
//  1
//  //!    *4: 1 for closed/isolated links
//  //!    f(F) for pipes
//  //!    f(F) for head pumps
//  //!    f(Hstart,Hend) for power pumps
//  //!    0 for active PRVs
//  //!    f(F) for open PRVs
//  //!    f(F) for open or active TCVs
//  //!    f(F) for open FCVs
//  //!    1 for active FCVs
//  //!    *5: 0 for inactive leaks
//  //!    0 for leaks at isolated junctions
//  //!    f(H-z) otherwise
//
//  void initialize_jacobian();
//};
}  // namespace linear_system
}  // namespace pipenetwork
#endif  // PIPE_NETWORK_MATRIX_ASSEMBLER_H
