#include "matrix_assembler.h"

// pipenetwork::MatrixAssembler::MatrixAssembler(
//        const std::shared_ptr<pipenetwork::Mesh>& mesh,
//        std::shared_ptr<pipenetwork::Curves>& curves_info, bool pdd_mode)
//        : mesh_{mesh}, curves_info_{curves_info}, pdd_{pdd_mode} {
//
//    init_variable_vectors();
////    assemble_balance_headloss_matrix();
////    initialize_jacobian();
//
//}

// Initialize variable vector
// 0 to nnode-1 element: nodal head
// nnode to 2*nnode-1 element: nodal demand
// 2*nnode to 2*nnode+nlinks-1 element: pipe discharge
// 2*nnode+nlinks to 2*nnode+nlinks+nleaks element: leak discharge
void pipenetwork::linear_system::Variables::init_variable_vectors() {
  auto nodes = mesh_->nodes();
  auto links = mesh_->links();
  auto leak_ids = mesh_->leak_nids();
  Index nnodes = nodes->nnodes();
  Index nlinks = links->nlinks();
  Index nleaks = leak_ids.size();

  resize_init_variables(nnodes, nlinks, nleaks);

  init_nodes_vecs(nodes);
  init_link_vecs(links, 2 * nnodes);
  init_leak_nodes(nodes, leak_ids, 2 * nnodes + nlinks);
}

void pipenetwork::linear_system::Variables::resize_init_variables(
    Index nnodes, Index nlinks, Index nleaks) {

  variable_vec_.resize(2 * nnodes + nlinks + nleaks);
  demands_heads_vec_.resize(nnodes);
  elevations_.resize(nnodes);
  link_resistance_coeff_vec_.resize(nlinks);
  link_minor_loss_coeff_vec_.resize(nlinks);
  leak_areas_.resize(nleaks);
}

void pipenetwork::linear_system::Variables::init_nodes_vecs(
    const std::shared_ptr<pipenetwork::MeshNodes>& nodes) {
  auto junctions_map = nodes->junctions();
  auto reservors_map = nodes->reservoirs();
  auto nnodes = nodes->nnodes();

  for (const auto& index_junction : junctions_map) {
    auto nid = index_junction.first;
    auto junction = index_junction.second;
    auto junc_prop = junction->property();
    demands_heads_vec_[nid] = junc_prop.demand;
    variable_vec_[nnodes + nid] = junc_prop.demand;

    elevations_[nid] = junc_prop.elevation;
    variable_vec_[nid] = junc_prop.elevation;
  }

  for (const auto& index_res : reservors_map) {
    auto nid = index_res.first;
    auto res = index_res.second;
    demands_heads_vec_[nid] = res->head();
    variable_vec_[nid] = res->head();
    elevations_[nid] = res->head();
    variable_vec_[nnodes + nid] = 0;
  }
}

void pipenetwork::linear_system::Variables::init_link_vecs(
    const std::shared_ptr<pipenetwork::MeshLinks>& links,
    Index link_start_idx) {
  auto pipes_map = links->pipes();
  auto pumps_map = links->pumps();
  auto valves_map = links->valves();
  auto links_map = links->links();

  update_minor_resis_coeff(pipes_map);
  update_minor_resis_coeff(pumps_map);
  update_minor_resis_coeff(valves_map);

  for (const auto& index_link : links_map) {
    auto lid = index_link.first;
    variable_vec_[link_start_idx + lid] = INIT_FLOWRATE;
  }
}

template <typename LinkMap>
void pipenetwork::linear_system::Variables::update_minor_resis_coeff(
    const LinkMap& linkmap) {

  for (const auto& index_link : linkmap) {
    auto lid = index_link.first;
    auto link = index_link.second;
    auto link_res_coeff = get_link_res_coeff(link);
    auto minor_loss_coeff = get_link_minor_coeff(link);

    link_resistance_coeff_vec_[lid] = link_res_coeff;
    link_minor_loss_coeff_vec_[lid] = minor_loss_coeff;
  }
}

double pipenetwork::linear_system::Variables::get_link_res_coeff(
    const std::shared_ptr<pipenetwork::Pipe>& pipe) {
  auto pipe_property = pipe->property();
  auto link_res_coeff = HW_COEFF * std::pow(pipe_property.roughness, -1.852) *
                        std::pow(pipe_property.diameter, -4.871) *
                        pipe_property.length;
  return link_res_coeff;
}

double pipenetwork::linear_system::Variables::get_link_res_coeff(
    const std::shared_ptr<pipenetwork::Valve>& valve) {

  auto valve_property = valve->property();
  double link_res_coeff = 0;
  if (valve_property.type == ValveType::TCVALVE) {
    link_res_coeff =
        8 * valve_property.setting /
        (G * std::pow(PI, 2) * std::pow(valve_property.diameter, 4));
  }
  return link_res_coeff;
}

void pipenetwork::linear_system::Variables::init_leak_nodes(
    const std::shared_ptr<pipenetwork::MeshNodes>& nodes,
    const std::vector<Index>& leak_ids, Index leak_start_idx) {
  auto junction_map = nodes->junctions();
  Index i = 0;
  for (const auto& leak_id : leak_ids) {
    auto leak_junction = junction_map.at(leak_id);
    leak_areas_[i] = leak_junction->leak_area();
    variable_vec_[leak_start_idx + i] = 0;
    i++;
  }
}

pipenetwork::linear_system::Residuals::Residuals(
    const std::shared_ptr<pipenetwork::Mesh>& mesh,
    const std::shared_ptr<pipenetwork::linear_system::Variables>& vars,
    const std::shared_ptr<Curves>& curves)
    : mesh_{mesh},
      vars_{vars},
      variable_vec_{vars_->variables_vec()},
      curves_info_{curves} {

  auto var_size = vars_->variables_vec().size();
  residual_vec_.resize(var_size);

  nnodes_ = mesh_->nodes()->nnodes();
  nlinks_ = mesh_->links()->nlinks();

  assemble_iso_masks();
  assemble_balance_headloss_matrix();
}

void pipenetwork::linear_system::Residuals::assemble_iso_masks() {
  iso_junctions_mask_.setZero(nnodes_);
  for (auto& nid : mesh_->iso_nodes()) {
    iso_junctions_mask_[nid] = 1;
  }

  iso_links_mask_.setZero(nlinks_);
  for (auto& lid : mesh_->iso_links()) {
    iso_links_mask_[lid] = 1;
  }

  Eigen::VectorXd ones_link;
  Eigen::VectorXd ones_nodes;
  ones_link.setOnes(nlinks_);
  ones_nodes.setOnes(nnodes_);
  connect_junctions_mask_ = ones_nodes - iso_junctions_mask_;
  connect_links_mask_ = ones_link - iso_links_mask_;
}

// Assemble node balance and link headloss matrix, which contains relation
// between discharges and heads. These matrix are parts of big jacobian matrix
// and can be used for fast residual computation
void pipenetwork::linear_system::Residuals::assemble_balance_headloss_matrix() {
  std::vector<Eigen::Triplet<double>> update_balance, update_headloss;
  auto links = mesh_->links();

  update_balance.reserve(nnodes_ + nlinks_);
  update_headloss.reserve(nnodes_ + nlinks_);

  for (const auto& index_link : links->links()) {
    auto lid = index_link.first;
    auto link = index_link.second;
    auto out_node_id = link->nodes().first->id();
    auto in_node_id = link->nodes().second->id();

    update_balance.emplace_back(out_node_id, lid, -1);
    update_balance.emplace_back(in_node_id, lid, 1);

    update_headloss.emplace_back(lid, out_node_id, 1);
    update_headloss.emplace_back(lid, in_node_id, -1);
  }

  node_balance_mat_.resize(nnodes_, nlinks_);
  node_balance_mat_.setFromTriplets(update_balance.begin(),
                                    update_balance.end());

  headloss_mat_.resize(nlinks_, nnodes_);
  headloss_mat_.setFromTriplets(update_headloss.begin(), update_headloss.end());
}

void pipenetwork::linear_system::Residuals::assemble_residual() {

  assemble_node_balance_residual();
  assemble_demand_head_residual();

  assemble_headloss_residual_pipe();
  assemble_headloss_residual_pump();
  assemble_headloss_residual_valve();

  assemble_leak_residual();
}

void pipenetwork::linear_system::Residuals::assemble_node_balance_residual() {
  residual_vec_.segment(0, nnodes_) =
      node_balance_mat_ * (variable_vec_.segment(2 * nnodes_, nlinks_)) -
      variable_vec_.segment(nnodes_, nnodes_);

  // correct on leaking nodes
  auto leak_ids = mesh_->leak_nids();
  int leak_idx = 0;
  for (const auto& leak_id : leak_ids) {
    residual_vec_.coeffRef(leak_id) -=
        variable_vec_.coeffRef(2 * nnodes_ + nlinks_ + leak_idx);
    ++leak_idx;
  }
}

void pipenetwork::linear_system::Residuals::assemble_demand_head_residual() {

  auto demands_heads_vec = vars_->demands_heads_vec();

  // demand for junctions
  residual_vec_.segment(nnodes_, nnodes_) =
      iso_junctions_mask_.array() * variable_vec_.segment(0, nnodes_).array() +
      connect_junctions_mask_.array() *
          (variable_vec_.segment(nnodes_, nnodes_) - demands_heads_vec).array();

  // correct residuals for sources (head for reservoir/tanks)
  auto nsrcs = mesh_->nodes()->nreservoirs();
  auto njunctions = mesh_->nodes()->njunctions();
  residual_vec_.segment(nnodes_ + njunctions, nsrcs) =
      variable_vec_.segment(njunctions, nsrcs) -
      demands_heads_vec.segment(njunctions, nsrcs);
}

void pipenetwork::linear_system::Residuals::assemble_headloss_residual_pipe() {

  auto npipes = mesh_->links()->npipes();
  auto sign_array = (variable_vec_.segment(2 * nnodes_, npipes))
                        .unaryExpr([](double x) {
                          //                            return (x < 0) ? -1 : 1;
                          if (x > 0) return 1.0;
                          return -1.0;
                        })
                        .array();  // get the sign of discharges
  // case 1, discharges that exceed the boundary of HW_Q2, use normal
  // hazen-william equation
  auto case1_bool = (variable_vec_.segment(2 * nnodes_, npipes))
                        .unaryExpr([](double x) {
                          if (std::abs(x) > HW_Q2) return 1.0;
                          return 0.0;
                        })
                        .array();
  // case 2, discharges that fall in between HW_Q1 and HW_Q2, use polynomial
  // approximation
  auto case2_bool = (variable_vec_.segment(2 * nnodes_, npipes))
                        .unaryExpr([](double x) {
                          if ((std::abs(x) < HW_Q2) && (std::abs(x) > HW_Q1))
                            return 1.0;
                          return 0.0;
                        })
                        .array();
  // case 3, discharges that are smaller than HW_Q1 , approximate 0 headloss in
  // this case
  auto case3_bool = (variable_vec_.segment(2 * nnodes_, npipes))
                        .unaryExpr([](double x) {
                          if (std::abs(x) < HW_Q1) return 1.0;
                          return 0.0;
                        })
                        .array();
  auto discharge_abs_array =
      ((variable_vec_.segment(2 * nnodes_, npipes)).array()).abs();
  auto head_diff_array = (headloss_mat_ * (variable_vec_.segment(0, nnodes_)))
                             .segment(0, npipes)
                             .array();

  auto hw_poly_vec = curves_info_->poly_coeffs()["HW_POLY_VEC"];

  auto link_resistance_coeff_vec = vars_->link_resistance_coeff_vec();

  residual_vec_.segment(2 * nnodes_, npipes) =
      (iso_links_mask_.segment(0, npipes).array() *
       variable_vec_.segment(2 * nnodes_, npipes).array()) +
      connect_links_mask_.segment(0, npipes).array() *
          (case1_bool *
               (sign_array *
                    link_resistance_coeff_vec.segment(0, npipes).array() *
                    discharge_abs_array.pow(1.852) -
                head_diff_array)

           + case2_bool *
                 (sign_array *
                      link_resistance_coeff_vec.segment(0, npipes).array() *
                      (hw_poly_vec[0] * discharge_abs_array.pow(3) +
                       hw_poly_vec[1] * discharge_abs_array.pow(2) +
                       hw_poly_vec[2] * discharge_abs_array + hw_poly_vec[3]) -
                  head_diff_array) +
           case3_bool *
               (sign_array *
                    link_resistance_coeff_vec.segment(0, npipes).array() *
                    HW_M * discharge_abs_array -
                head_diff_array));
}

void pipenetwork::linear_system::Residuals::assemble_headloss_residual_pump() {
  auto pumps_map = mesh_->links()->pumps();
  for (const auto& index_pump : pumps_map) {
    auto lid = index_pump.first;
    auto pump = index_pump.second;
    auto nodes = pump->nodes();
    auto pump_info = pump->property();
    auto link_flow = variable_vec_[2 * nnodes_ + lid];

    double pump_residual = 0;
    if (iso_links_mask_[lid] == 1) {
      pump_residual = link_flow;
    } else if (pump_info.type == PumpType::HEADPUMP) {
      auto pump_headgain = get_pump_headgain(pump, link_flow);
      pump_residual = pump_headgain - (variable_vec_[nodes.first->id()] -
                                       variable_vec_[nodes.second->id()]);

    } else if (pump_info.type == PumpType::POWERPUMP) {
      auto head_diff_array =
          (headloss_mat_ * (variable_vec_.segment(0, nnodes_))).array();
      pump_residual =
          pump_info.power + (head_diff_array[lid]) * link_flow * G * 1000.0;
    }

    residual_vec_[2 * nnodes_ + lid] = pump_residual;
  }
}

double pipenetwork::linear_system::Residuals::get_pump_headgain(
    const std::shared_ptr<pipenetwork::Pump>& pump, double link_flow) {
  auto pump_info = pump->property();

  auto pump_curve = curves_info_->pump_curves().at(
      curves_info_->pump_int_str(pump_info.curve_id));
  auto curve_coeff = pump_curve.head_curve_coefficients;

  double pump_headgain = 0;

  if (curve_coeff[2] > 1) {
    auto line_coeff = pump_curve.line_param;
    if (link_flow >= line_coeff[0]) {
      pump_headgain =
          curve_coeff[0] - curve_coeff[1] * std::pow(link_flow, curve_coeff[2]);
    } else {
      pump_headgain = PUMP_M * (link_flow - line_coeff[0]) + line_coeff[1];
    }
  } else {
    if (link_flow < PUMP_Q1) {
      pump_headgain = PUMP_M * link_flow + curve_coeff[0];
    } else if (link_flow < PUMP_Q2) {
      auto curve_poly_coeff = pump_curve.poly_coefficients;
      pump_headgain = curve_poly_coeff[0] * std::pow(link_flow, 3) +
                      curve_poly_coeff[1] * std::pow(link_flow, 2) +
                      curve_poly_coeff[2] * link_flow + curve_poly_coeff[3];
    } else {
      pump_headgain =
          curve_coeff[0] - curve_coeff[1] * std::pow(link_flow, curve_coeff[2]);
    }
  }
  return pump_headgain;
}

void pipenetwork::linear_system::Residuals::assemble_headloss_residual_valve() {

  auto valve_map = mesh_->links()->valves();
  for (const auto& index_valve : valve_map) {
    auto lid = index_valve.first;
    auto valve = index_valve.second;

    auto valve_info = valve->property();
    auto link_flow = variable_vec_[2 * nnodes_ + lid];

    double valve_residual = 0;
    if (iso_links_mask_[lid] == 1 || valve_info.status == LinkStatus::CLOSED) {
      valve_residual = link_flow;
    } else if (valve_info.status == LinkStatus::OPEN) {
      valve_residual = get_open_valve_residual(valve, link_flow);
    }
    if (valve_info.status == LinkStatus::ACTIVE) {
      valve_residual = get_active_valve_residual(valve, link_flow);
    }

    residual_vec_[2 * nnodes_ + lid] = valve_residual;
  }
}

double pipenetwork::linear_system::Residuals::get_open_valve_residual(
    const std::shared_ptr<Valve>& valve, double link_flow) {

  auto lid = valve->id();
  auto link_minor_loss_coeff_vec = vars_->link_minor_loss_coeff_vec();
  auto pipe_headloss = link_minor_loss_coeff_vec[lid] * std::pow(link_flow, 2);
  if (link_flow < 0) {
    pipe_headloss = -1 * pipe_headloss;
  }
  auto head_diff_array =
      (headloss_mat_ * (variable_vec_.segment(0, nnodes_))).array();

  double valve_residual = pipe_headloss - head_diff_array[lid];
  return valve_residual;
}

double pipenetwork::linear_system::Residuals::get_active_valve_residual(
    const std::shared_ptr<Valve>& valve, double link_flow) {
  auto valve_info = valve->property();
  auto nodes = valve->nodes();

  double valve_res = 0;
  if (valve_info.type == ValveType::PRVALVE) {
    auto end_node_idx = nodes.second->id();
    auto end_node_elevation = vars_->elevations()[end_node_idx];
    valve_res =
        variable_vec_[end_node_idx] - (valve_info.setting + end_node_elevation);
  } else if (valve_info.type == ValveType::FCVALVE) {
    valve_res = link_flow - valve_info.setting;
  } else if (valve_info.type == ValveType::TCVALVE) {
    auto link_resistance_coeff_vec = vars_->link_resistance_coeff_vec();
    auto coeff = link_resistance_coeff_vec[valve->id()];
    auto pipe_headloss = coeff * std::pow(link_flow, 2);
    if (link_flow < 0) {
      pipe_headloss = -1 * pipe_headloss;
    }
    auto head_diff_array =
        (headloss_mat_ * (variable_vec_.segment(0, nnodes_))).array();

    valve_res = pipe_headloss - head_diff_array[valve->id()];
  }
  return valve_res;
}

// Residual for leak from holes. piece-wise function, use for loop... (number of
// leak nodes is usually small too)

void pipenetwork::linear_system::Residuals::assemble_leak_residual() {
  double m = 1e-4;
  Index leak_idx = 2 * nnodes_ + nlinks_;
  auto leak_ids = mesh_->leak_nids();
  auto elevations = vars_->elevations();
  auto leak_areas = vars_->leak_areas();

  double leak_residual = 0;
  for (const auto& leak_id : leak_ids) {
    if (iso_junctions_mask_[leak_id] == 1) {
      leak_residual = variable_vec_[leak_idx];
    } else {
      auto p = variable_vec_[leak_id] - elevations[leak_id];
      // case 1, no pressure no leak
      if (p < m) {
        leak_residual = variable_vec_[leak_idx] - m * p;
      }
      // case 3, normal leak equation
      else {
        auto i = leak_idx - 2 * nnodes_ - nlinks_;
        leak_residual = variable_vec_[leak_idx] -
                        LEAK_COEFF * leak_areas[i] * std::sqrt(2 * G * p);
      }
      residual_vec_[leak_idx] = leak_residual;
      ++leak_idx;
    }
  }
}
