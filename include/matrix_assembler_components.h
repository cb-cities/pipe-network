#ifndef PIPE_NETWORK_MATRIX_ASSEMBLER_COMPONENTS_H
#define PIPE_NETWORK_MATRIX_ASSEMBLER_COMPONENTS_H

#include "curves.h"
#include "mesh.h"
#include "settings.h"

namespace pipenetwork {
namespace linear_system {

//! Vairables class
//! \brief Class for linear system variables
class Variables {
 public:
  explicit Variables(const std::shared_ptr<Mesh>& mesh) : mesh_{mesh} {
    assemble_iso_masks();
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

  //! Get the isolated junctions vector (boolean mask purpose)
  const Eigen::VectorXd& iso_junctions_mask() const {
    return iso_junctions_mask_;
  }

  //! Get the connected junctions vector (boolean mask purpose)
  const Eigen::VectorXd& connect_junctions_mask() const {
    return connect_junctions_mask_;
  }

  //! Get the isolated/closed links vector (boolean mask purpose)
  const Eigen::VectorXd& iso_links_mask() const { return iso_links_mask_; }

  //! Get the connected links vector (boolean mask purpose)
  const Eigen::VectorXd& connect_links_mask() const {
    return connect_links_mask_;
  }

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

struct HW_vectors {
  Eigen::ArrayXd sign_array;
  Eigen::ArrayXd case1_bool;
  Eigen::ArrayXd case2_bool;
  Eigen::ArrayXd case3_bool;
  Eigen::ArrayXd discharge_abs_array;
  Eigen::ArrayXd head_diff_array;
  Eigen::Vector4d hw_poly_vec;
};
struct PDD_vectors {
  Eigen::ArrayXd pressure;
  Eigen::ArrayXd case1_bool;
  Eigen::ArrayXd case2_bool;
  Eigen::ArrayXd case3_bool;
  Eigen::ArrayXd case4_bool;
  Eigen::ArrayXd case5_bool;
  Eigen::Vector4d pdd1_poly_vec;
  Eigen::Vector4d pdd2_poly_vec;
};

//! Residuals class
//! \brief Class for linear system residuals
class Residuals {
 public:
  Residuals(const std::shared_ptr<Mesh>& mesh,
            const std::shared_ptr<Variables>& vars,
            const std::shared_ptr<Curves>& curves);

  //! assemble the residual vector for demand driven simulation
  void assemble_residual();

  //! assemble the residual vector for pressure demand driven simulation (PDD)
  void assemble_residual_pdd();

  //! get the residual vector
  const Eigen::VectorXd& residual_vec() const { return residual_vec_; }

  //! get the harzan-williams vectors
  const HW_vectors& hw_vectors() const { return hw_vec_; }

  //! get the pressure demand driven vectors
  const PDD_vectors& pdd_vectors() const { return pdd_vec_; }

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
  //! Variable vector ptr
  Eigen::VectorXd& variable_vec_;
  //! number of nodes
  Index nnodes_;
  //! number of links
  Index nlinks_;

  //! information for Hazan-williams equation
  HW_vectors hw_vec_;
  //! information for pressure demand driven (PDD) demand head equation
  PDD_vectors pdd_vec_;
  //! assemble the balance headloss matrix for fast residual computation
  void assemble_balance_headloss_matrix();

  //! assemble node balance residual (0-nnodes)
  void assemble_node_balance_residual();

  //! assemble demand head residual (nnodes-2*nnodes)
  void assemble_demand_head_residual();

  //! assemble demand head residual for pressure demand driven mode
  void assemble_demand_head_residual_pdd();

  //! assemble necessary information for pressure demand equation (PDD)
  void update_pdd_vectors();

  //! assemble headloss residual for pipes (hazan-williams equation)
  void assemble_headloss_residual_pipe();
  //! assemble necessary information for headloss (hazan-williams equation)
  void update_hw_vectors();

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

//! Jacobian class
//! \brief Class for linear system jacobian
class Jacobian {
 public:
  Jacobian(const std::shared_ptr<Mesh>& mesh,
           const std::shared_ptr<Variables>& vars,
           const std::shared_ptr<Residuals>& residuals,
           const std::shared_ptr<Curves>& curves);

  //! Method to update jacobian matrix from the variable vector (Demand driven
  //! mode)
  void update_jacobian();
  //! Method to update jacobian matrix from the variable vector (Pressure demand
  //! driven mode)
  void update_jacobian_pdd();

  //! Method to get jacobian matrix
  const Eigen::SparseMatrix<double, Eigen::RowMajor>& jac_matrix() const {
    return jac_;
  }

 private:
  //! the mesh ptr
  std::shared_ptr<Mesh> mesh_;
  //! curves ptr
  std::shared_ptr<Residuals> residuals_;
  //! curves ptr
  std::shared_ptr<Curves> curves_info_;
  //! the vriable ptr
  std::shared_ptr<Variables> vars_;
  //! Variable vector ptr
  Eigen::VectorXd& variable_vec_;
  //! Jacobian matrix
  Eigen::SparseMatrix<double, Eigen::RowMajor> jac_;
  //! map of for sub-jacobian triplets (row, col, value)
  std::map<std::string, std::vector<Eigen::Triplet<double>>> sub_jac_trip_;
  //! number of nodes
  Index nnodes_;
  //! number of links
  Index nlinks_;
  //! Initialize the jacobian matrix, the Jacobian matrix is splitted into 9
  //! submatrices: sub jacobian A: node_bal equation with respect to demand sub
  //! jacobian B: node_bal equation with respect to flow sub jacobian C:
  //! node_bal equation with respect to  leak flow sub jacobian D: demand/head
  //! equation with respect to head sub jacobian E: demand/head equation with
  //! respect to flow sub jacobian F: headloss equation with respect to head sub
  //! jacobian G: headloss equation with respect to flow sub jacobian H: leak
  //! flow to flow sub jacobian I: leak flow to leak flow
  void initialize_jacobian();
  void initialize_subjacA();
  void initialize_subjacBnF();
  void initialize_subjacC();
  void initialize_subjacDnE();
  void initialize_subjacG();
  void initialize_subjacHnI();

  //! set jacobian constant based on connectivity
  void set_jac_const();
  void set_jacF_const();

  // Jac_d: pressure-demand equation
  void update_jac_d();
  // Jac_f: for power pump only
  void update_jac_f();
  // jac_g: headloss equation (harzen-william)
  void update_jac_g_pipe();
  void update_jac_g_pump();
  void update_jac_g_valve();
  // helper function to get pump jacobian
  double get_pump_jac(const std::shared_ptr<pipenetwork::Pump>& pump,
                      double link_flow);
  // jac_h: leak equation
  void update_jac_h();
};
}  // namespace linear_system
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_MATRIX_ASSEMBLER_COMPONENTS_H
