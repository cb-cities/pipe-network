#include "matrix_assembler_components.h"

void pipenetwork::linear_system::Variables::assemble_iso_masks() {
  auto nnodes = mesh_->nodes()->nnodes();
  auto nlinks = mesh_->links()->nlinks();

  iso_junctions_mask_.setZero(nnodes);
  for (auto& nid : mesh_->iso_nodes()) {
    iso_junctions_mask_[nid] = 1;
  }

  iso_links_mask_.setZero(nlinks);
  for (auto& lid : mesh_->iso_links()) {
    iso_links_mask_[lid] = 1;
  }

  Eigen::VectorXd ones_link;
  Eigen::VectorXd ones_nodes;
  ones_link.setOnes(nlinks);
  ones_nodes.setOnes(nnodes);
  connect_junctions_mask_ = ones_nodes - iso_junctions_mask_;
  connect_links_mask_ = ones_link - iso_links_mask_;
}

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

  assemble_balance_headloss_matrix();
  assemble_residual_pdd();  // initialize HW/PDD vectors for jacobian matrix
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

void pipenetwork::linear_system::Residuals::assemble_residual_pdd() {

  assemble_residual();
  assemble_demand_head_residual_pdd();
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
  auto iso_junction_mask = vars_->iso_junctions_mask();
  auto connect_junctions_mask = vars_->connect_junctions_mask();
  // demand for junctions
  residual_vec_.segment(nnodes_, nnodes_) =
      iso_junction_mask.array() * variable_vec_.segment(0, nnodes_).array() +
      connect_junctions_mask.array() *
          (variable_vec_.segment(nnodes_, nnodes_) - demands_heads_vec).array();

  // correct residuals for sources (head for reservoir/tanks)
  auto nsrcs = mesh_->nodes()->nreservoirs();
  auto njunctions = mesh_->nodes()->njunctions();
  residual_vec_.segment(nnodes_ + njunctions, nsrcs) =
      variable_vec_.segment(njunctions, nsrcs) -
      demands_heads_vec.segment(njunctions, nsrcs);
}

void pipenetwork::linear_system::Residuals::update_pdd_vectors() {
  auto pressure = (variable_vec_.segment(0, nnodes_) - vars_->elevations());
  pdd_vec_.pressure = pressure;

  // case 1, pressure smaller than min pressure, no water
  pdd_vec_.case1_bool = pressure
                            .unaryExpr([](double x) {
                              if (x <= MIN_PRESSURE) return 1.0;
                              return 0.0;
                            })
                            .array();
  // case 2, pressure larger than min pressure but in a small range, use
  // polynomial approximation
  pdd_vec_.case2_bool =
      pressure
          .unaryExpr([](double x) {
            if ((x > MIN_PRESSURE) && (x < (MIN_PRESSURE + PDD_DELTA)))
              return 1.0;
            return 0.0;
          })
          .array();
  // case 3, pressure close to normal pressure, use polynomial approximation
  pdd_vec_.case3_bool =
      pressure
          .unaryExpr([](double x) {
            if ((x > (NORMAL_PRESSURE - PDD_DELTA)) && (x < (NORMAL_PRESSURE)))
              return 1.0;
            return 0.0;
          })
          .array();
  // case 4, pressure above normal pressure, demand can be met
  pdd_vec_.case4_bool = pressure
                            .unaryExpr([](double x) {
                              if ((x > NORMAL_PRESSURE)) return 1.0;
                              return 0.0;
                            })
                            .array();
  // case 5, pressure falls in between min pressure and normal pressure, use
  // pressure-demand equation
  pdd_vec_.case5_bool = pressure
                            .unaryExpr([](double x) {
                              if ((x > (MIN_PRESSURE + PDD_DELTA)) &&
                                  (x < (NORMAL_PRESSURE - PDD_DELTA)))
                                return 1.0;
                              return 0.0;
                            })
                            .array();
  pdd_vec_.pdd1_poly_vec = curves_info_->poly_coeffs()["PDD_POLY_VEC1"];
  pdd_vec_.pdd2_poly_vec = curves_info_->poly_coeffs()["PDD_POLY_VEC2"];
}

void pipenetwork::linear_system::Residuals::
    assemble_demand_head_residual_pdd() {

  update_pdd_vectors();

  auto iso_links_mask = vars_->iso_links_mask();
  auto connect_links_mask = vars_->connect_links_mask();
  auto demands_heads_vec = vars_->demands_heads_vec();

  residual_vec_.segment(nnodes_, nnodes_) =
      iso_links_mask.array() * variable_vec_.segment(0, nnodes_).array() +
      connect_links_mask.array() *
          (pdd_vec_.case1_bool *
               (variable_vec_.segment(nnodes_, nnodes_).array() -
                demands_heads_vec.array() * PDD_SLOPE *
                    (pdd_vec_.pressure.array() - MIN_PRESSURE)) +
           pdd_vec_.case2_bool *
               ((variable_vec_.segment(nnodes_, nnodes_).array()) -
                demands_heads_vec.array() *
                    (pdd_vec_.pdd1_poly_vec[0] *
                         pdd_vec_.pressure.array().pow(3) +
                     pdd_vec_.pdd1_poly_vec[1] *
                         pdd_vec_.pressure.array().pow(2) +
                     pdd_vec_.pdd1_poly_vec[2] * pdd_vec_.pressure.array() +
                     pdd_vec_.pdd1_poly_vec[3])) +
           pdd_vec_.case3_bool *
               ((variable_vec_.segment(nnodes_, nnodes_).array()) -
                demands_heads_vec.array() *
                    (pdd_vec_.pdd2_poly_vec[0] *
                         pdd_vec_.pressure.array().pow(3) +
                     pdd_vec_.pdd2_poly_vec[1] *
                         pdd_vec_.pressure.array().pow(2) +
                     pdd_vec_.pdd2_poly_vec[2] * pdd_vec_.pressure.array() +
                     pdd_vec_.pdd2_poly_vec[3])) +
           pdd_vec_.case4_bool *
               ((variable_vec_.segment(nnodes_, nnodes_).array()) -
                demands_heads_vec.array() *
                    (PDD_SLOPE * (pdd_vec_.pressure.array() - NORMAL_PRESSURE) +
                     1)) +
           pdd_vec_.case5_bool *
               ((variable_vec_.segment(nnodes_, nnodes_).array()) -
                demands_heads_vec.array() *
                    ((pdd_vec_.pressure.array().abs() - MIN_PRESSURE) /
                     (NORMAL_PRESSURE - MIN_PRESSURE))
                        .abs()
                        .pow(0.5)));

  // correct residuals for sources (head for reservoir/tanks)
  auto nsrcs = mesh_->nodes()->nreservoirs();
  auto njunctions = mesh_->nodes()->njunctions();
  residual_vec_.segment(nnodes_ + njunctions, nsrcs) =
      variable_vec_.segment(njunctions, nsrcs) -
      demands_heads_vec.segment(njunctions, nsrcs);
}

void pipenetwork::linear_system::Residuals::update_hw_vectors() {
  auto npipes = mesh_->links()->npipes();
  hw_vec_.sign_array = (variable_vec_.segment(2 * nnodes_, npipes))
                           .unaryExpr([](double x) {
                             //                            return (x < 0) ? -1 :
                             //                            1;
                             if (x > 0) return 1.0;
                             return -1.0;
                           })
                           .array();  // get the sign of discharges
  // case 1, discharges that exceed the boundary of HW_Q2, use normal
  // hazen-william equation
  hw_vec_.case1_bool = (variable_vec_.segment(2 * nnodes_, npipes))
                           .unaryExpr([](double x) {
                             if (std::abs(x) > HW_Q2) return 1.0;
                             return 0.0;
                           })
                           .array();
  // case 2, discharges that fall in between HW_Q1 and HW_Q2, use polynomial
  // approximation
  hw_vec_.case2_bool = (variable_vec_.segment(2 * nnodes_, npipes))
                           .unaryExpr([](double x) {
                             if ((std::abs(x) < HW_Q2) && (std::abs(x) > HW_Q1))
                               return 1.0;
                             return 0.0;
                           })
                           .array();
  // case 3, discharges that are smaller than HW_Q1 , approximate 0 headloss in
  // this case
  hw_vec_.case3_bool = (variable_vec_.segment(2 * nnodes_, npipes))
                           .unaryExpr([](double x) {
                             if (std::abs(x) < HW_Q1) return 1.0;
                             return 0.0;
                           })
                           .array();
  hw_vec_.discharge_abs_array =
      ((variable_vec_.segment(2 * nnodes_, npipes)).array()).abs();
  hw_vec_.head_diff_array =
      (headloss_mat_ * (variable_vec_.segment(0, nnodes_)))
          .segment(0, npipes)
          .array();
  hw_vec_.hw_poly_vec = curves_info_->poly_coeffs()["HW_POLY_VEC"];
}
void pipenetwork::linear_system::Residuals::assemble_headloss_residual_pipe() {

  auto npipes = mesh_->links()->npipes();
  auto link_resistance_coeff_vec = vars_->link_resistance_coeff_vec();
  auto iso_links_mask = vars_->iso_links_mask();
  auto connect_links_mask = vars_->connect_links_mask();

  update_hw_vectors();

  residual_vec_.segment(2 * nnodes_, npipes) =
      (iso_links_mask.segment(0, npipes).array() *
       variable_vec_.segment(2 * nnodes_, npipes).array()) +
      connect_links_mask.segment(0, npipes).array() *
          (hw_vec_.case1_bool *
               (hw_vec_.sign_array *
                    link_resistance_coeff_vec.segment(0, npipes).array() *
                    hw_vec_.discharge_abs_array.pow(1.852) -
                hw_vec_.head_diff_array)

           + hw_vec_.case2_bool *
                 (hw_vec_.sign_array *
                      link_resistance_coeff_vec.segment(0, npipes).array() *
                      (hw_vec_.hw_poly_vec[0] *
                           hw_vec_.discharge_abs_array.pow(3) +
                       hw_vec_.hw_poly_vec[1] *
                           hw_vec_.discharge_abs_array.pow(2) +
                       hw_vec_.hw_poly_vec[2] * hw_vec_.discharge_abs_array +
                       hw_vec_.hw_poly_vec[3]) -
                  hw_vec_.head_diff_array) +
           hw_vec_.case3_bool *
               (hw_vec_.sign_array *
                    link_resistance_coeff_vec.segment(0, npipes).array() *
                    HW_M * hw_vec_.discharge_abs_array -
                hw_vec_.head_diff_array));
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
    if (vars_->iso_links_mask()[lid] == 1) {
      pump_residual = link_flow;
    } else if (pump_info.type == PumpType::HEADPUMP) {
      auto pump_headgain = get_pump_headgain(pump, link_flow);
      pump_residual = pump_headgain - (variable_vec_[nodes.second->id()] -
                                       variable_vec_[nodes.first->id()]);

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
    if (vars_->iso_links_mask()[lid] == 1 ||
        valve_info.status == LinkStatus::CLOSED) {
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
    if (vars_->iso_junctions_mask()[leak_id] == 1) {
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

pipenetwork::linear_system::Jacobian::Jacobian(
    const std::shared_ptr<pipenetwork::Mesh>& mesh,
    const std::shared_ptr<pipenetwork::linear_system::Variables>& vars,
    const std::shared_ptr<pipenetwork::linear_system::Residuals>& residuals,
    const std::shared_ptr<Curves>& curves)
    : mesh_{mesh},
      vars_{vars},
      residuals_{residuals},
      curves_info_{curves},
      variable_vec_{vars_->variables_vec()} {

  nnodes_ = mesh_->nodes()->nnodes();
  nlinks_ = mesh_->links()->nlinks();

  initialize_jacobian();
  update_jacobian();
}

void pipenetwork::linear_system::Jacobian::initialize_jacobian() {
  initialize_subjacA();
  initialize_subjacBnF();
  initialize_subjacC();
  initialize_subjacDnE();
  initialize_subjacG();
  initialize_subjacHnI();

  // loop through all the sub-matrix to form the initial jacobian matrix
  auto leaks_id = mesh_->leak_nids();
  jac_.resize(2 * nnodes_ + nlinks_ + leaks_id.size(),
              2 * nnodes_ + nlinks_ + leaks_id.size());

  std::vector<Eigen::Triplet<double>> jac_triplet;
  jac_triplet.reserve(10 * nnodes_);
  for (auto const& trip_map : sub_jac_trip_) {
    jac_triplet.insert(std::end(jac_triplet), std::begin(trip_map.second),
                       std::end(trip_map.second));
  }
  jac_.setFromTriplets(jac_triplet.begin(), jac_triplet.end());
}

void pipenetwork::linear_system::Jacobian::initialize_subjacA() {
  std::vector<Eigen::Triplet<double>> update;
  update.reserve(nnodes_);
  for (int row_idx = 0; row_idx < nnodes_; row_idx++) {
    auto col_idx = row_idx + nnodes_;
    update.emplace_back(row_idx, col_idx, -1);
  }
  sub_jac_trip_["jac_a"] = update;
}

void pipenetwork::linear_system::Jacobian::initialize_subjacBnF() {
  std::vector<Eigen::Triplet<double>> update_jacb, update_jacf;
  update_jacb.reserve(2 * nlinks_);
  update_jacf.reserve(2 * nlinks_);

  auto links = mesh_->links();
  auto pump_start_lid = links->npipes();
  auto pump_end_lid = pump_start_lid + links->npumps();
  for (const auto& index_link : links->links()) {
    auto lid = index_link.first;
    auto link = index_link.second;
    auto out_node_id = link->nodes().first->id();
    auto in_node_id = link->nodes().second->id();

    update_jacb.emplace_back(out_node_id, lid + 2 * nnodes_, -1);
    update_jacb.emplace_back(in_node_id, lid + 2 * nnodes_, 1);

    double valout = -1, valin = 1;
    // reverse on pumps
    if (lid >= pump_start_lid && lid < pump_end_lid) {
      valin *= -1;
      valout *= -1;
    }
    update_jacf.emplace_back(lid + 2 * nnodes_, out_node_id, valout);
    update_jacf.emplace_back(lid + 2 * nnodes_, in_node_id, valin);
  }

  sub_jac_trip_["jac_b"] = update_jacb;
  sub_jac_trip_["jac_f"] = update_jacf;
}

void pipenetwork::linear_system::Jacobian::initialize_subjacC() {
  auto leak_ids = mesh_->leak_nids();

  std::vector<Eigen::Triplet<double>> update;
  update.reserve(leak_ids.size());

  Index idx = 0;
  for (const auto& leak_id : leak_ids) {
    auto col_idx = 2 * nnodes_ + nlinks_ + idx;
    update.emplace_back(leak_id, col_idx, -1);
    ++idx;
  }
  sub_jac_trip_["jac_c"] = update;
}

void pipenetwork::linear_system::Jacobian::initialize_subjacDnE() {
  auto nodes = mesh_->nodes();
  auto src_start_nid = nodes->njunctions();
  auto src_end_nid = src_start_nid + nodes->nreservoirs();

  std::vector<Eigen::Triplet<double>> updateD, updateE;
  updateD.reserve(nnodes_);
  updateE.reserve(nnodes_);

  for (int row_idx = nnodes_; row_idx < nnodes_ + nnodes_; row_idx++) {
    auto nid = row_idx - nnodes_;
    if (nid >= src_start_nid && nid < src_end_nid) {
      updateD.emplace_back(row_idx, row_idx - nnodes_, 1);  // sources
      updateE.emplace_back(row_idx, row_idx, 0);
    } else {
      updateD.emplace_back(row_idx, row_idx - nnodes_, 0);
      updateE.emplace_back(row_idx, row_idx, 1);
    }
  }
  sub_jac_trip_["jac_d"] = updateD;
  sub_jac_trip_["jac_e"] = updateE;
}

void pipenetwork::linear_system::Jacobian::initialize_subjacG() {
  std::vector<Eigen::Triplet<double>> update;
  update.reserve(nlinks_);

  for (int row_idx = 0; row_idx < nlinks_; row_idx++) {
    auto col_idx = row_idx;
    update.emplace_back(row_idx + 2 * nnodes_, col_idx + 2 * nnodes_, 1);
  }
  sub_jac_trip_["jac_g"] = update;
}

void pipenetwork::linear_system::Jacobian::initialize_subjacHnI() {
  auto leak_ids = mesh_->leak_nids();

  std::vector<Eigen::Triplet<double>> updateH, updateI;
  updateH.reserve(leak_ids.size());
  updateI.reserve(leak_ids.size());

  unsigned idx = 0;
  for (const auto& leak_id : leak_ids) {
    auto col_idx = leak_id;
    auto row_idx = 2 * nnodes_ + nlinks_ + idx;
    updateH.emplace_back(row_idx, col_idx, 0);
    updateI.emplace_back(row_idx, row_idx, 1);
    ++idx;
  }
  sub_jac_trip_["jac_h"] = updateH;
  sub_jac_trip_["jac_i"] = updateI;
}

void pipenetwork::linear_system::Jacobian::update_jacobian() {
  set_jac_const();
  // Jac_f: for power pump only
  update_jac_f();
  // jac_g: headloss equation (harzen-william)
  update_jac_g_pipe();
  update_jac_g_pump();
  update_jac_g_valve();
  // jac_h: leak equation
  update_jac_h();
}

void pipenetwork::linear_system::Jacobian::update_jacobian_pdd() {
  update_jacobian();
  // Jac_d: pressure-demand equation
  update_jac_d();
}

void pipenetwork::linear_system::Jacobian::set_jac_const() {
  // get trip f and isolated junctions+sources

  Eigen::VectorXd iso_junc_src = vars_->iso_junctions_mask();
  Eigen::VectorXd connect_junc_no_src = vars_->connect_junctions_mask();

  auto reservoir_map = mesh_->nodes()->reservoirs();
  for (const auto& index_src : reservoir_map) {
    auto nid = index_src.first;
    iso_junc_src[nid] = 1;
    connect_junc_no_src[nid] = 0;
  }

  // change jacobian matrix entries based on connection status
  int count = 0;
  for (const auto& trip_d : sub_jac_trip_["jac_d"]) {
    jac_.coeffRef(trip_d.row(), trip_d.col()) = iso_junc_src[count];
    ++count;
  }
  count = 0;
  for (const auto& trip_e : sub_jac_trip_["jac_e"]) {
    jac_.coeffRef(trip_e.row(), trip_e.col()) = connect_junc_no_src[count];
    ++count;
  }

  set_jacF_const();
}

void pipenetwork::linear_system::Jacobian::set_jacF_const() {
  auto trip_f = sub_jac_trip_["jac_f"];
  auto connect_links = vars_->connect_links_mask();
  for (int i = 0; i < nlinks_; ++i) {
    jac_.coeffRef(trip_f[2 * i].row(), trip_f[2 * i].col()) =
        connect_links[i] * trip_f[2 * i].value();
    jac_.coeffRef(trip_f[2 * i + 1].row(), trip_f[2 * i + 1].col()) =
        connect_links[i] * trip_f[2 * i + 1].value();
  }
  // iterate through all the valves for jac_f corrections
  auto valves_map = mesh_->links()->valves();
  auto npipes = mesh_->links()->npipes();
  auto npumps = mesh_->links()->npumps();

  for (const auto& index_valve : valves_map) {
    Index valve_idx = index_valve.first;
    auto valve = index_valve.second;
    auto valve_info = valve->property();

    if (valve_info.status == LinkStatus::ACTIVE) {
      if (valve_info.type == ValveType::PRVALVE) {
        jac_.coeffRef(trip_f[2 * valve_idx].row(),
                      trip_f[2 * valve_idx].col()) = 0;
      } else if (valve_info.type == ValveType::FCVALVE) {
        jac_.coeffRef(trip_f[2 * valve_idx].row(),
                      trip_f[2 * valve_idx].col()) = 0;
        jac_.coeffRef(trip_f[2 * valve_idx + 1].row(),
                      trip_f[2 * valve_idx + 1].col()) = 0;
      }
    }
  }
}
//
void pipenetwork::linear_system::Jacobian::update_jac_d() {
  auto pdd_vec = residuals_->pdd_vectors();

  auto iso_junctions = vars_->iso_junctions_mask();
  auto connect_junctions = vars_->connect_junctions_mask();
  auto demans_head_vec = vars_->demands_heads_vec();

  auto vals =
      iso_junctions.array() +
      connect_junctions.array() *
          (pdd_vec.case1_bool * (-PDD_SLOPE * demans_head_vec.array() *
                                 variable_vec_.segment(0, nnodes_).array()) +
           pdd_vec.case2_bool * (-demans_head_vec.array() *
                                 (3 * pdd_vec.pdd1_poly_vec[0] *
                                      pdd_vec.pressure.array().pow(2) +
                                  2 * pdd_vec.pdd1_poly_vec[1] *
                                      pdd_vec.pressure.array().pow(1) +
                                  pdd_vec.pdd1_poly_vec[2])) +
           pdd_vec.case3_bool * (-demans_head_vec.array() *
                                 (3 * pdd_vec.pdd2_poly_vec[0] *
                                      pdd_vec.pressure.array().pow(2) +
                                  2 * pdd_vec.pdd2_poly_vec[1] *
                                      pdd_vec.pressure.array().pow(1) +
                                  pdd_vec.pdd2_poly_vec[2])) +
           pdd_vec.case4_bool * (-PDD_SLOPE * demans_head_vec.array() *
                                 variable_vec_.segment(0, nnodes_).array()) +
           pdd_vec.case5_bool *
               (-0.5 * demans_head_vec.array() /
                (NORMAL_PRESSURE - MIN_PRESSURE) *
                ((pdd_vec.pressure.array().abs() - MIN_PRESSURE) /
                 (NORMAL_PRESSURE - MIN_PRESSURE))
                    .abs()
                    .pow(-0.5)));

  auto trip_d = sub_jac_trip_["jac_d"];
  auto njunctions = mesh_->nodes()->njunctions();
  for (Index i = 0; i < njunctions; ++i) {
    jac_.coeffRef(trip_d[i].row(), trip_d[i].col()) = vals[i];
  }
}

void pipenetwork::linear_system::Jacobian::update_jac_f() {
  auto trip_f = sub_jac_trip_["jac_f"];
  auto pumps_map = mesh_->links()->pumps();
  auto iso_links = vars_->iso_links_mask();

  for (const auto& index_pump : pumps_map) {
    auto lid = index_pump.first;
    auto pump = index_pump.second;

    if (iso_links[lid] == 0) {
      auto pump_info = pump->property();

      if (pump_info.type == PumpType::POWERPUMP) {
        auto link_flow = variable_vec_[2 * nnodes_ + lid];

        jac_.coeffRef(trip_f[2 * lid].row(), trip_f[2 * lid].col()) =
            1000.0 * G * link_flow;
        jac_.coeffRef(trip_f[2 * lid + 1].row(), trip_f[2 * lid + 1].col()) =
            -1000.0 * G * link_flow;
      }
    }
  }
}

void pipenetwork::linear_system::Jacobian::update_jac_g_pipe() {
  auto npipes = mesh_->links()->npipes();
  auto link_resistance_coeff_vec = vars_->link_resistance_coeff_vec();
  auto iso_links_mask = vars_->iso_links_mask();
  auto connect_links_mask = vars_->connect_links_mask();
  auto hw_vectors = residuals_->hw_vectors();

  auto vals =
      (iso_links_mask.segment(0, npipes)).array() +
      connect_links_mask.segment(0, npipes).array() *
          (hw_vectors.case1_bool * 1.852 *
               link_resistance_coeff_vec.segment(0, npipes).array() *
               hw_vectors.discharge_abs_array.pow(.852) +
           hw_vectors.case2_bool *
               link_resistance_coeff_vec.segment(0, npipes).array() *
               (3 * hw_vectors.hw_poly_vec[0] *
                    hw_vectors.discharge_abs_array.pow(2) +
                2 * hw_vectors.hw_poly_vec[1] * hw_vectors.discharge_abs_array +
                1 * hw_vectors.hw_poly_vec[2]) +
           hw_vectors.case3_bool *
               link_resistance_coeff_vec.segment(0, npipes).array() * HW_M);
  auto trip_g = sub_jac_trip_["jac_g"];
  for (int i = 0; i < npipes; ++i) {
    jac_.coeffRef(trip_g[i].row(), trip_g[i].col()) = vals[i];
  }
}

void pipenetwork::linear_system::Jacobian::update_jac_g_pump() {

  auto trip_g = sub_jac_trip_["jac_g"];
  auto pumps_map = mesh_->links()->pumps();

  for (const auto& index_pump : pumps_map) {
    auto lid = index_pump.first;
    auto pump = index_pump.second;
    auto nodes = pump->nodes();
    auto pump_info = pump->property();
    auto link_flow = variable_vec_[2 * nnodes_ + lid];

    double pump_jac = 0;
    if (vars_->iso_links_mask()[lid] == 1) {
      pump_jac = 1;
    } else if (pump_info.type == PumpType::HEADPUMP) {
      pump_jac = get_pump_jac(pump, link_flow);

    } else if (pump_info.type == PumpType::POWERPUMP) {
      auto start_node_idx = nodes.first->id();
      auto end_node_idx = nodes.second->id();
      pump_jac = (1000.0 * G * variable_vec_[start_node_idx] -
                  1000.0 * G * variable_vec_[end_node_idx]);
    }

    jac_.coeffRef(trip_g[lid].row(), trip_g[lid].col()) = pump_jac;
  }
}

double pipenetwork::linear_system::Jacobian::get_pump_jac(
    const std::shared_ptr<pipenetwork::Pump>& pump, double link_flow) {
  auto pump_info = pump->property();

  auto pump_curve = curves_info_->pump_curves().at(
      curves_info_->pump_int_str(pump_info.curve_id));
  auto curve_coeff = pump_curve.head_curve_coefficients;

  double pump_jac = 0;

  if (curve_coeff[2] > 1) {
    auto line_coeff = pump_curve.line_param;
    if (link_flow >= line_coeff[0]) {
      pump_jac = -curve_coeff[1] * curve_coeff[2] *
                 std::pow(link_flow, (curve_coeff[2] - 1));

    } else {
      pump_jac = PUMP_M;
    }
  } else {
    if (link_flow < PUMP_Q1) {
      pump_jac = PUMP_M;
    } else if (link_flow < PUMP_Q2) {
      auto curve_poly_coeff = pump_curve.poly_coefficients;
      pump_jac = 3 * curve_poly_coeff[0] * std::pow(link_flow, 2) +
                 2 * curve_poly_coeff[1] * link_flow + curve_poly_coeff[2];
    } else {
      pump_jac = -curve_coeff[1] * curve_coeff[2] *
                 std::pow(link_flow, (curve_coeff[2] - 1));
    }
  }
  return pump_jac;
}

void pipenetwork::linear_system::Jacobian::update_jac_g_valve() {

  auto hw_vectors = residuals_->hw_vectors();
  auto link_resistance_coeff_vec = vars_->link_resistance_coeff_vec();
  auto link_minor_loss_coeff_vec = vars_->link_minor_loss_coeff_vec();

  auto head_diff_array = hw_vectors.head_diff_array;
  auto trip_g = sub_jac_trip_["jac_g"];

  // iterate through all the valves for jac_f corrections
  auto valves_map = mesh_->links()->valves();
  auto npipes = mesh_->links()->npipes();
  auto npumps = mesh_->links()->npumps();
  auto iso_links = vars_->iso_links_mask();

  for (const auto& index_valve : valves_map) {
    Index lid = index_valve.first;
    auto valve = index_valve.second;

    auto valve_info = valve->property();
    if (iso_links[lid] == 1) {
      jac_.coeffRef(trip_g[lid].row(), trip_g[lid].col()) = 1;
    } else {
      double val = 0;
      auto link_flow = variable_vec_[2 * nnodes_ + lid];

      // PRV
      if (valve_info.type == ValveType::PRVALVE) {
        if (valve_info.status == LinkStatus::CLOSED) {
          val = 1;
        } else if (valve_info.status == LinkStatus::ACTIVE) {
          val = 0;
        }
      }

      // FCV
      else if (valve_info.type == ValveType::FCVALVE) {
        if (valve_info.status == LinkStatus::CLOSED) {
          val = 1;
        } else if (valve_info.status == LinkStatus::ACTIVE) {
          val = 1;
        }
      }

      // TCV
      else if (valve_info.type == ValveType::TCVALVE) {
        if (valve_info.status == LinkStatus::CLOSED) {
          val = link_flow;
        } else if (valve_info.status == LinkStatus::ACTIVE) {
          auto coeff = link_resistance_coeff_vec[lid];
          val = 2 * coeff * std::abs(link_flow);
        }
      }

      if (valve_info.status == LinkStatus::OPEN) {
        auto coeff = link_minor_loss_coeff_vec[lid];
        val = 2 * coeff * std::abs(link_flow);
      }
      jac_.coeffRef(trip_g[lid].row(), trip_g[lid].col()) = val;
    }
  }
}

void pipenetwork::linear_system::Jacobian::update_jac_h() {
  double m = 1e-4;
  Index leak_idx = 2 * nnodes_ + nlinks_;
  auto leak_ids = mesh_->leak_nids();
  auto elevations = vars_->elevations();
  auto leak_areas = vars_->leak_areas();
  auto trip_h = sub_jac_trip_["jac_h"];

  double val = 0;
  for (const auto& leak_id : leak_ids) {
    auto p = variable_vec_[leak_id] - elevations[leak_id];
    // case 1, no pressure no leak
    if (p < m) {
      val = -m;
    }
    // case 3, normal leak equation
    else {
      auto i = leak_idx - 2 * nnodes_ - nlinks_;
      val = -0.5 * LEAK_COEFF * leak_areas[i] * std::sqrt(2 * G) *
            std::pow(p, -0.5);
    }
    jac_.coeffRef(trip_h[leak_idx - 2 * nnodes_ - nlinks_].row(),
                  trip_h[leak_idx - 2 * nnodes_ - nlinks_].col()) = val;
    ++leak_idx;
  }
}
