#include "matrix_assembler.h"

// Initialize variable vector
// 1 to nnode element: nodal head
// nnode+1 to 2*nnode element: nodal demand
// 2*nnode+1 to 2*nnode+npipe element: pipe discharge
// 2*nnode+npipe+1 to 3*nnode+npipe element: leak discharge
void pipenetwork::MatrixAssembler::init_variable_vector() {
  variable_vec_->resize(3 * nnodes_ + nlinks_);
  demands_heads_vec_.resize(nnodes_);
  elevations_.resize(nnodes_);
  link_resistance_coeff_vec_.resize(nlinks_);
  link_minor_loss_coeff_vec_.resize(nlinks_);

  Index idx = 0, leak_idx = 0;
  for (const auto& node : mesh_->nodes()) {
    Index index_nd = idx + nnodes_;
    Index index_nl = leak_idx + 2 * nnodes_ + nlinks_;

    // get node information, determine possible leak nodes and assemble demands
    // heads vector
    auto info = node.second->nodal_info();
    switch (static_cast<int>(info["type"])) {
      case JUNCTION:
        demands_heads_vec_[idx] = info["demand"];
        elevations_[idx] = info["elevation"];
        // leak information
        if (info["leak_area"] > 0) {
          std::string leak_node_name = node.second->id();
          double leak_area = info["leak_area"];
          leak_ids_.emplace_back(leak_node_name);
          leak_area_.emplace_back(leak_area);
          // register leak pdd curve
          curves_info_->add_leak_poly_vec(leak_node_name, leak_area);
          // nodal leak discharge vector
          variable_vec_->coeffRef(index_nl) = node.second->sim_leak();
          ++leak_idx;
        }
        break;
      default:
        demands_heads_vec_[idx] = info["head"];
        source_idx_.emplace_back(idx);
    }

    // nodal head vector
    variable_vec_->coeffRef(idx) = node.second->sim_head();
    // nodal demand vector
    variable_vec_->coeffRef(index_nd) = node.second->sim_demand();
    // id map
    node_id_map_.emplace(node.second->id(), idx);
    node_idx_map_.emplace(idx, node.second->id());

    ++idx;
  }

  // loop through links
  for (const auto& link : mesh_->links()) {
    Index index_pd = idx + nnodes_;
    // discharge vector
    variable_vec_->coeffRef(index_pd) = link->sim_discharge();
    // id map
    link_id_map_.emplace(link->id(), idx - nnodes_);
    link_idx_map_.emplace(idx - nnodes_, link->id());

    auto info = link->link_info();
    double val = 0;
    double minor_loss_coeff = 0;
    // construct the energy loss coefficients for pipes and valves
    if (info["type"] != POWERPUMP && info["type"] != HEADPUMP) {
      minor_loss_coeff = 8 * info["minor_loss"] /
                         (G * std::pow(PI, 2) * std::pow(info["diameter"], 4));
      switch (static_cast<pipenetwork::Link_type>(info["type"])) {
        // link resistence coefficients for pipes (Harzan-williams)
        case PIPE:
          val = HW_COEFF * std::pow(info["roughness"], -1.852) *
                std::pow(info["diameter"], -4.871) * info["length"];
          break;
        case TCVALVE:
          val = 8 * info["setting"] /
                (G * std::pow(PI, 2) * std::pow(info["diameter"], 4));

          break;
      }
    }

    link_minor_loss_coeff_vec_[idx - nnodes_] = minor_loss_coeff;
    link_resistance_coeff_vec_[idx - nnodes_] = val;
    ++idx;
  }
  // use only part of the original array due to uncertainty of number of leak
  // nodes
  variable_vec_ = std::make_shared<Eigen::VectorXd>(
      variable_vec_->segment(0, leak_idx + 2 * nnodes_ + nlinks_));
}

// Assemble node balance and link headloss matrix, which contains relation
// between discharges and heads. These matrix are parts of big jacobian matrix
// and can be used for fast residual computation
void pipenetwork::MatrixAssembler::assemble_balance_headloss_matrix() {

  std::vector<Eigen::Triplet<double>> update_balance, update_headloss,
      update_jacb, update_jacf;
  update_balance.reserve(nnodes_ + nlinks_);
  update_headloss.reserve(nnodes_ + nlinks_);
  update_jacb.reserve(nnodes_ + nlinks_);
  update_jacf.reserve(nnodes_ + nlinks_);

  for (const auto& link : mesh_->links()) {
    auto out_node_id = link->nodes().first->id();
    auto in_node_id = link->nodes().second->id();
    auto link_type = link->link_info()["type"];

    update_balance.emplace_back(node_id_map_[out_node_id],
                                link_id_map_[link->id()], -1);
    update_balance.emplace_back(node_id_map_[in_node_id],
                                link_id_map_[link->id()], 1);

    update_jacb.emplace_back(node_id_map_[out_node_id],
                             link_id_map_[link->id()] + 2 * nnodes_, -1);
    update_jacb.emplace_back(node_id_map_[in_node_id],
                             link_id_map_[link->id()] + 2 * nnodes_, 1);

    update_headloss.emplace_back(link_id_map_[link->id()],
                                 node_id_map_[out_node_id], 1);
    update_headloss.emplace_back(link_id_map_[link->id()],
                                 node_id_map_[in_node_id], -1);

    double valout = -1, valin = 1;

    if (link_type == HEADPUMP || link_type == POWERPUMP) {
      valout = 1;
      valin = -1;
    }

    update_jacf.emplace_back(link_id_map_[link->id()] + 2 * nnodes_,
                             node_id_map_[out_node_id], valout);
    update_jacf.emplace_back(link_id_map_[link->id()] + 2 * nnodes_,
                             node_id_map_[in_node_id], valin);
  }

  node_balance_mat_.resize(nnodes_, nlinks_);
  node_balance_mat_.setFromTriplets(update_balance.begin(),
                                    update_balance.end());

  headloss_mat_.resize(nlinks_, nnodes_);
  headloss_mat_.setFromTriplets(update_headloss.begin(), update_headloss.end());

  sub_jac_trip_["jac_b"] = update_jacb;
  sub_jac_trip_["jac_f"] = update_jacf;

  //  std::cout << (*node_balance_mat_) << std::endl;
}

// Calculate and assemble residual vector
// Nodal balance equation:
// Residual = flow = Inflow - Outflow - Demand
//
// demand_pressure equation :
// For junctions: Residual = demand_var-demand_desired (demand driven mode)
// Residual = demand_var - pressure_demand_func(pressure demand driven mode)
// For sources: Residual = head_var-head_desired
//
// Headloss equation (Hazen-Williams):
// Residual = \frac{10.67 \times length
//            \times pipe_discharge^1.852}{pipe_roughness^1.852 \times
//            (2radius)^4.8704}- (Start_node_head - end_node_head)

// The network has nnode Nodal balance equation ,nnode demand_pressure equation,
// nnode node leak equation and npipe Headloss equation Thus the residual vector
// has (3*nnode+npipe) elements
void pipenetwork::MatrixAssembler::assemble_residual() {
  residual_vec_->resize(variable_vec_->size());

  // node balance residual
  residual_vec_->segment(0, nnodes_) =
      node_balance_mat_ * (variable_vec_->segment(2 * nnodes_, nlinks_)) -
      variable_vec_->segment(nnodes_, nnodes_);

  // get leak effect
  int leak_count = 0;
  for (const auto& leak_id : leak_ids_) {
    auto idx = node_id_map_[leak_id];
    residual_vec_->coeffRef(idx) -=
        variable_vec_->coeffRef(2 * nnodes_ + nlinks_ + leak_count);
    ++leak_count;
  }

  // demand or head residual
  assemble_demand_head_residual();

  // headloss residual
  assemble_headloss_residual_pipe();
  assemble_headloss_residual_pump();
  assemble_headloss_residual_valve();
  // leak residual
  assemble_leak_residual();
}

// Residual for mass conservation equation (what comes in must go out)
void pipenetwork::MatrixAssembler::assemble_demand_head_residual() {
  if (!pdd_) {
    // demand for junctions
    residual_vec_->segment(nnodes_, nnodes_) =
        iso_junctions_.array() * variable_vec_->segment(0, nnodes_).array() +
        connect_junctions_.array() *
            (variable_vec_->segment(nnodes_, nnodes_) - demands_heads_vec_)
                .array();

  } else {
    auto pressure = (variable_vec_->segment(0, nnodes_) - elevations_);

    // case 1, pressure smaller than min pressure, no water
    auto case1_bool = pressure
                          .unaryExpr([](double x) {
                            if (x <= MIN_PRESSURE) return 1.0;
                            return 0.0;
                          })
                          .array();
    // case 2, pressure larger than min pressure but in a small range, use
    // polynomial approximation
    auto case2_bool =
        pressure
            .unaryExpr([](double x) {
              if ((x > MIN_PRESSURE) && (x < (MIN_PRESSURE + PDD_DELTA)))
                return 1.0;
              return 0.0;
            })
            .array();
    // case 3, pressure close to normal pressure, use polynomial approximation
    auto case3_bool = pressure
                          .unaryExpr([](double x) {
                            if ((x > (NORMAL_PRESSURE - PDD_DELTA)) &&
                                (x < (NORMAL_PRESSURE)))
                              return 1.0;
                            return 0.0;
                          })
                          .array();
    // case 4, pressure above normal pressure, demand can be met
    auto case4_bool = pressure
                          .unaryExpr([](double x) {
                            if ((x > NORMAL_PRESSURE)) return 1.0;
                            return 0.0;
                          })
                          .array();
    // case 5, pressure falls in between min pressure and normal pressure, use
    // pressure-demand equation
    auto case5_bool = pressure
                          .unaryExpr([](double x) {
                            if ((x > (MIN_PRESSURE + PDD_DELTA)) &&
                                (x < (NORMAL_PRESSURE - PDD_DELTA)))
                              return 1.0;
                            return 0.0;
                          })
                          .array();
    auto pdd1_poly_vec = curves_info_->poly_coeffs()["PDD_POLY_VEC1"];
    auto pdd2_poly_vec = curves_info_->poly_coeffs()["PDD_POLY_VEC2"];

    residual_vec_->segment(nnodes_, nnodes_) =
        iso_junctions_.array() * variable_vec_->segment(0, nnodes_).array() +
        connect_junctions_.array() *
            (case1_bool * (variable_vec_->segment(nnodes_, nnodes_).array() -
                           demands_heads_vec_.array() * PDD_SLOPE *
                               (pressure.array() - MIN_PRESSURE)) +
             case2_bool * ((variable_vec_->segment(nnodes_, nnodes_).array()) -
                           demands_heads_vec_.array() *
                               (pdd1_poly_vec[0] * pressure.array().pow(3) +
                                pdd1_poly_vec[1] * pressure.array().pow(2) +
                                pdd1_poly_vec[2] * pressure.array() +
                                pdd1_poly_vec[3])) +
             case3_bool * ((variable_vec_->segment(nnodes_, nnodes_).array()) -
                           demands_heads_vec_.array() *
                               (pdd2_poly_vec[0] * pressure.array().pow(3) +
                                pdd2_poly_vec[1] * pressure.array().pow(2) +
                                pdd2_poly_vec[2] * pressure.array() +
                                pdd2_poly_vec[3])) +
             case4_bool *
                 ((variable_vec_->segment(nnodes_, nnodes_).array()) -
                  demands_heads_vec_.array() *
                      (PDD_SLOPE * (pressure.array() - NORMAL_PRESSURE) + 1)) +
             case5_bool * ((variable_vec_->segment(nnodes_, nnodes_).array()) -
                           demands_heads_vec_.array() *
                               ((pressure.array().abs() - MIN_PRESSURE) /
                                (NORMAL_PRESSURE - MIN_PRESSURE))
                                   .abs()
                                   .pow(0.5)));
  }
  // correct residuals for sources (head for reservoir/tanks)
  for (const auto& idx : source_idx_) {
    (*residual_vec_)[idx + nnodes_] =
        (*variable_vec_)[idx] - demands_heads_vec_[idx];
  }
  //  std::cout<<"demand_residual" << residual_vec_->segment(nnodes_,
  //  nnodes_).norm()<<std::endl;
}

// Residual for leak from holes. piece-wise function, use for loop... (number of
// leak nodes is usually small too)
void pipenetwork::MatrixAssembler::assemble_leak_residual() {
  double m = 1e-11;
  Index leak_idx = 2 * nnodes_ + nlinks_;
  for (const auto& leak_id : leak_ids_) {
    auto i = leak_idx - 2 * nnodes_ - nlinks_;
    auto idx = node_id_map_[leak_id];
    if (iso_junctions_[idx] == 1) {
      residual_vec_->coeffRef(leak_idx) = (*variable_vec_)[leak_idx];
    } else {
      auto p = (*variable_vec_)[idx] - elevations_[idx];

      // case 1, no pressure no leak
      if (p < m) {
        residual_vec_->coeffRef(leak_idx) = (*variable_vec_)[leak_idx] - m * p;
      }
      // case 2, around the boundary, use polynomial approximation
      else if (p < 1e-4) {
        auto leak_poly_coef = curves_info_->poly_coeffs()[leak_id];
        residual_vec_->coeffRef(leak_idx) =
            (*variable_vec_)[leak_idx] - leak_poly_coef[0] * std::pow(p, 3) -
            leak_poly_coef[1] * std::pow(p, 2) - leak_poly_coef[2] * p -
            leak_poly_coef[3];
      }
      // case 3, normal leak equation
      else {
        residual_vec_->coeffRef(leak_idx) =
            (*variable_vec_)[leak_idx] -
            LEAK_COEFF * leak_area_[i] * std::sqrt(2 * G * p);
      }
      ++leak_idx;
    }
  }
}

// residual for energy conservation equation (Hazen-williams). Piece-wise
// function for stability concern
void pipenetwork::MatrixAssembler::assemble_headloss_residual_pipe() {

  auto sign_array = (variable_vec_->segment(2 * nnodes_, npipes_))
                        .unaryExpr([](double x) {
                          //                            return (x < 0) ? -1 : 1;
                          if (x > 0) return 1.0;
                          return -1.0;
                        })
                        .array();  // get the sign of discharges
  // case 1, discharges that exceed the boundary of HW_Q2, use normal
  // hazen-william equation
  auto case1_bool = (variable_vec_->segment(2 * nnodes_, npipes_))
                        .unaryExpr([](double x) {
                          if (std::abs(x) > HW_Q2) return 1.0;
                          return 0.0;
                        })
                        .array();
  // case 2, discharges that fall in between HW_Q1 and HW_Q2, use polynomial
  // approximation
  auto case2_bool = (variable_vec_->segment(2 * nnodes_, npipes_))
                        .unaryExpr([](double x) {
                          if ((std::abs(x) < HW_Q2) && (std::abs(x) > HW_Q1))
                            return 1.0;
                          return 0.0;
                        })
                        .array();
  // case 3, discharges that are smaller than HW_Q1 , approximate 0 headloss in
  // this case
  auto case3_bool = (variable_vec_->segment(2 * nnodes_, npipes_))
                        .unaryExpr([](double x) {
                          if (std::abs(x) < HW_Q1) return 1.0;
                          return 0.0;
                        })
                        .array();
  auto discharge_abs_array =
      ((variable_vec_->segment(2 * nnodes_, npipes_)).array()).abs();
  auto head_diff_array = (headloss_mat_ * (variable_vec_->segment(0, nnodes_)))
                             .segment(0, npipes_)
                             .array();

  auto hw_poly_vec = curves_info_->poly_coeffs()["HW_POLY_VEC"];

  residual_vec_->segment(2 * nnodes_, npipes_) =
      (iso_links_.segment(0, npipes_).array() *
       variable_vec_->segment(2 * nnodes_, npipes_).array()) +
      connect_links_.segment(0, npipes_).array() *
          (case1_bool *
               (sign_array *
                    link_resistance_coeff_vec_.segment(0, npipes_).array() *
                    discharge_abs_array.pow(1.852) -
                head_diff_array)

           + case2_bool *
                 (sign_array *
                      link_resistance_coeff_vec_.segment(0, npipes_).array() *
                      (hw_poly_vec[0] * discharge_abs_array.pow(3) +
                       hw_poly_vec[1] * discharge_abs_array.pow(2) +
                       hw_poly_vec[2] * discharge_abs_array + hw_poly_vec[3]) -
                  head_diff_array) +
           case3_bool *
               (sign_array *
                    link_resistance_coeff_vec_.segment(0, npipes_).array() *
                    HW_M * discharge_abs_array -
                head_diff_array));
}

void pipenetwork::MatrixAssembler::assemble_headloss_residual_pump() {
  // iterate through all the pumps
  for (int i = 0; i < npumps_; ++i) {
    int pump_idx = npipes_ + i;
    auto link_flow = (*variable_vec_)[2 * nnodes_ + pump_idx];
    if (iso_links_[pump_idx] == 1) {
      residual_vec_->array()[2 * nnodes_ + pump_idx] = link_flow;
    } else {
      auto pump = mesh_->links().at(pump_idx);
      auto pump_info = pump->link_info();
      auto nodes = pump->nodes();
      auto start_node_idx = node_id_map_.at(nodes.first->id());
      auto end_node_idx = node_id_map_.at(nodes.second->id());

      // head pump
      if (pump_info["type"] == HEADPUMP) {
        double pump_headgain;
        auto pump_curve = curves_info_->pump_curves().at(
            curves_info_->pump_int_str(pump_info.at("curve_name")));
        auto curve_coeff = pump_curve.head_curve_coefficients;

        if (curve_coeff[2] > 1) {
          auto line_coeff = pump_curve.line_param;
          if (link_flow >= line_coeff[0]) {
            pump_headgain =
                curve_coeff[0] -
                curve_coeff[1] * std::pow(link_flow, curve_coeff[2]);

          }

          else {
            pump_headgain =
                PUMP_M * (link_flow - line_coeff[0]) + line_coeff[1];
          }
        } else {
          if (link_flow < PUMP_Q1) {
            pump_headgain = PUMP_M * link_flow + curve_coeff[0];
          } else if (link_flow < PUMP_Q2) {
            auto curve_poly_coeff = pump_curve.poly_coefficients;
            pump_headgain = curve_poly_coeff[0] * std::pow(link_flow, 3) +
                            curve_poly_coeff[1] * std::pow(link_flow, 2) +
                            curve_poly_coeff[2] * link_flow +
                            curve_poly_coeff[3];
          } else {
            pump_headgain =
                curve_coeff[0] -
                curve_coeff[1] * std::pow(link_flow, curve_coeff[2]);
          }
        }
        residual_vec_->array()[2 * nnodes_ + pump_idx] =
            pump_headgain - (variable_vec_->array()[end_node_idx] -
                             variable_vec_->array()[start_node_idx]);

      }
      // power pump
      else if (pump_info["type"] == POWERPUMP) {
        auto head_diff_array =
            (headloss_mat_ * (variable_vec_->segment(0, nnodes_))).array();

        residual_vec_->array()[2 * nnodes_ + pump_idx] =
            pump_info["power"] +
            (head_diff_array[pump_idx]) * link_flow * G * 1000.0;
      }
    }
  }
}

void pipenetwork::MatrixAssembler::assemble_headloss_residual_valve() {
  auto head_diff_array =
      (headloss_mat_ * (variable_vec_->segment(0, nnodes_))).array();

  // iterate through all the valves
  for (int i = 0; i < nvalves_; ++i) {
    int valve_idx = npipes_ + npumps_ + i;
    auto link_flow = (*variable_vec_)[2 * nnodes_ + valve_idx];
    if (iso_links_[valve_idx] == 1) {
      residual_vec_->array()[2 * nnodes_ + valve_idx] = link_flow;
    } else {

      auto valve = mesh_->links().at(valve_idx);
      auto valve_info = valve->link_info();
      auto nodes = valve->nodes();
      auto start_node_idx = node_id_map_.at(nodes.first->id());
      auto end_node_idx = node_id_map_.at(nodes.second->id());
      auto end_node_elevation = nodes.second->nodal_info()["elevation"];

      double val = 0;

      // PRV
      if (valve_info["type"] == PRVALVE) {
        if ((valve->link_status() == CLOSED)) {
          val = link_flow;
        } else if ((valve->link_status() == ACTIVE)) {
          val = variable_vec_->array()[end_node_idx] -
                (valve_info["setting"] + end_node_elevation);
        } else if ((valve->link_status() == OPEN)) {
          auto coeff = link_minor_loss_coeff_vec_[valve_idx];
          auto pipe_headloss = coeff * std::pow(link_flow, 2);
          if (link_flow < 0) {
            pipe_headloss = -1 * pipe_headloss;
          }
          val = pipe_headloss - head_diff_array[valve_idx];
        }

      }
      // FCV
      else if (valve_info["type"] == FCVALVE) {
        if ((valve->link_status() == CLOSED)) {
          val = link_flow;
        } else if ((valve->link_status() == ACTIVE)) {
          val = link_flow - valve_info["setting"];
        } else if ((valve->link_status() == OPEN)) {
          auto coeff = link_minor_loss_coeff_vec_[valve_idx];
          auto pipe_headloss = coeff * std::pow(link_flow, 2);
          if (link_flow < 0) {
            pipe_headloss = -1 * pipe_headloss;
          }
          val = pipe_headloss - head_diff_array[valve_idx];
        }
      }
      // TCV
      else if (valve_info["type"] == TCVALVE) {
        if ((valve->link_status() == CLOSED)) {
          val = link_flow;
        } else if ((valve->link_status() == ACTIVE)) {
          auto coeff = link_resistance_coeff_vec_[valve_idx];
          auto pipe_headloss = coeff * std::pow(link_flow, 2);
          if (link_flow < 0) {
            pipe_headloss = -1 * pipe_headloss;
          }
          val = pipe_headloss - head_diff_array[valve_idx];
        } else if ((valve->link_status() == OPEN)) {
          auto coeff = link_minor_loss_coeff_vec_[valve_idx];
          auto pipe_headloss = coeff * std::pow(link_flow, 2);
          if (link_flow < 0) {
            pipe_headloss = -1 * pipe_headloss;
          }
          val = pipe_headloss - head_diff_array[valve_idx];
        }
      }
      residual_vec_->array()[2 * nnodes_ + valve_idx] = val;
    }
  }
}

// initialize the jacobian matrix, and store sub-jacobian information for quick
// update during the iteration process
void pipenetwork::MatrixAssembler::initialize_jacobian() {
  std::vector<Eigen::Triplet<double>> update;
  update.reserve(nnodes_ + nlinks_);
  // JacA:
  for (int row_idx = 0; row_idx < nnodes_; row_idx++) {
    auto col_idx = row_idx + nnodes_;
    update.emplace_back(row_idx, col_idx, -1);
  }
  sub_jac_trip_["jac_a"] = update;
  update.clear();

  // JacB is set from the assemble_balance_headloss_matrix function

  // JacC
  Index idx = 0;
  for (const auto& leak_id : leak_ids_) {
    auto row_idx = node_id_map_[leak_id];
    auto col_idx = 2 * nnodes_ + nlinks_ + idx;
    update.emplace_back(row_idx, col_idx, -1);
    ++idx;
  }
  sub_jac_trip_["jac_c"] = update;
  update.clear();

  // JacD
  for (int row_idx = nnodes_; row_idx < nnodes_ + nnodes_; row_idx++) {
    auto col_idx = row_idx - nnodes_;
    if (std::find(source_idx_.begin(), source_idx_.end(), col_idx) !=
        source_idx_.end()) {
      update.emplace_back(row_idx, col_idx, 1);
    } else {
      update.emplace_back(row_idx, col_idx, 0);
    }
  }
  sub_jac_trip_["jac_d"] = update;
  update.clear();

  // JacE
  for (int row_idx = nnodes_; row_idx < nnodes_ + nnodes_; row_idx++) {
    auto col_idx = row_idx;
    if (std::find(source_idx_.begin(), source_idx_.end(), col_idx - nnodes_) !=
        source_idx_.end()) {
      update.emplace_back(row_idx, col_idx, 0);
    } else {
      update.emplace_back(row_idx, col_idx, 1);
    }
  }
  sub_jac_trip_["jac_e"] = update;
  update.clear();

  // JacF is set from the assemble_balance_headloss_matrix function

  // JacG
  for (int row_idx = 0; row_idx < nlinks_; row_idx++) {
    auto col_idx = row_idx;
    update.emplace_back(row_idx + 2 * nnodes_, col_idx + 2 * nnodes_, 1);
  }
  sub_jac_trip_["jac_g"] = update;
  update.clear();

  // JacH
  idx = 0;
  for (const auto& leak_id : leak_ids_) {
    auto col_idx = node_id_map_[leak_id];
    auto row_idx = 2 * nnodes_ + nlinks_ + idx;
    update.emplace_back(row_idx, col_idx, 0);
    ++idx;
  }
  sub_jac_trip_["jac_h"] = update;
  update.clear();

  // JacI
  idx = 0;
  for (const auto& leak_id : leak_ids_) {
    auto col_idx = 2 * nnodes_ + nlinks_ + idx;
    auto row_idx = 2 * nnodes_ + nlinks_ + idx;
    update.emplace_back(row_idx, col_idx, 1);
    ++idx;
  }
  sub_jac_trip_["jac_i"] = update;
  update.clear();

  // loop through all the sub-matrix to form the initial jacobian matrix

  jac_->resize(2 * nnodes_ + nlinks_ + leak_ids_.size(),
               2 * nnodes_ + nlinks_ + leak_ids_.size());
  std::vector<Eigen::Triplet<double>> jac_triplet;
  jac_triplet.reserve(10 * nnodes_);
  for (auto const& x : sub_jac_trip_) {
    //    std::cout << x.first  // string (key)
    //              << std::endl;
    jac_triplet.insert(std::end(jac_triplet), std::begin(x.second),
                       std::end(x.second));
  }
  jac_->setFromTriplets(jac_triplet.begin(), jac_triplet.end());
}

// update the jacobian matrix from updated variable vector
void pipenetwork::MatrixAssembler::update_jacobian() {
  set_jac_const();
  // Jac_d: pressure-demand equation
  update_jac_d();
  // Jac_f: for power pump only
  update_jac_f();
  // jac_g: headloss equation (harzen-william)
  update_jac_g_pipe();
  update_jac_g_pump();
  update_jac_g_valve();
  // jac_h: leak equation
  update_jac_h();
}

void pipenetwork::MatrixAssembler::set_jac_const() {
  // get trip f and isolated junctions+sources
  auto trip_f = sub_jac_trip_["jac_f"];
  auto iso_junc_src = iso_junctions_;
  auto connect_junc_no_src = connect_junctions_;
  for (const auto& src : source_idx_) {
    iso_junc_src[src] = 1;
    connect_junc_no_src[src] = 0;
  }

  // change jacobian matrix entries based on connection status
  int count = 0;
  if (!pdd_) {
    for (const auto& trip_d : sub_jac_trip_["jac_d"]) {
      jac_->coeffRef(trip_d.row(), trip_d.col()) = iso_junc_src[count];
      ++count;
    }
  }
  count = 0;
  for (const auto& trip_e : sub_jac_trip_["jac_e"]) {
    jac_->coeffRef(trip_e.row(), trip_e.col()) = connect_junc_no_src[count];
    ++count;
  }

  for (int i = 0; i < nlinks_; ++i) {
    jac_->coeffRef(trip_f[2 * i].row(), trip_f[2 * i].col()) =
        connect_links_[i] * trip_f[2 * i].value();
    jac_->coeffRef(trip_f[2 * i + 1].row(), trip_f[2 * i + 1].col()) =
        connect_links_[i] * trip_f[2 * i + 1].value();
  }

  // iterate through all the valves for jac_f corrections
  for (int i = 0; i < nvalves_; ++i) {
    int valve_idx = npipes_ + npumps_ + i;
    auto valve = mesh_->links().at(valve_idx);
    auto valve_info = valve->link_info();
    if ((valve->link_status() == ACTIVE)) {
      if (valve_info["type"] == PRVALVE) {
        jac_->coeffRef(trip_f[2 * valve_idx].row(),
                       trip_f[2 * valve_idx].col()) = 0;
      } else if (valve_info["type"] == FCVALVE) {
        jac_->coeffRef(trip_f[2 * valve_idx].row(),
                       trip_f[2 * valve_idx].col()) = 0;
        jac_->coeffRef(trip_f[2 * valve_idx + 1].row(),
                       trip_f[2 * valve_idx + 1].col()) = 0;
      }
    }
  }
}

// jac_d: Derivative of demand_pressure equation with respect to nodal head,
//            0 for junctions, 1 for reservoir/tank of demand driven model,
//            nonlinear for pressure_demand model.
void pipenetwork::MatrixAssembler::update_jac_d() {
  if (pdd_) {
    auto pressure = (variable_vec_->segment(0, nnodes_) - elevations_);
    // case 1, pressure smaller than min pressure, no water
    auto case1_bool = pressure
                          .unaryExpr([](double x) {
                            if (x <= MIN_PRESSURE) return 1.0;
                            return 0.0;
                          })
                          .array();
    // case 2, pressure larger than min pressure but in a small range, use
    // polynomial approximation
    auto case2_bool =
        pressure
            .unaryExpr([](double x) {
              if ((x > MIN_PRESSURE) && (x < (MIN_PRESSURE + PDD_DELTA)))
                return 1.0;
              return 0.0;
            })
            .array();
    // case 3, pressure close to normal pressure, use polynomial approximation
    auto case3_bool = pressure
                          .unaryExpr([](double x) {
                            if ((x > (NORMAL_PRESSURE - PDD_DELTA)) &&
                                (x < (NORMAL_PRESSURE)))
                              return 1.0;
                            return 0.0;
                          })
                          .array();
    // case 4, pressure above normal pressure, demand can be met
    auto case4_bool = pressure
                          .unaryExpr([](double x) {
                            if ((x > NORMAL_PRESSURE)) return 1.0;
                            return 0.0;
                          })
                          .array();
    // case 5, pressure falls in between min pressure and normal pressure, use
    // pressure-demand equation
    auto case5_bool = pressure
                          .unaryExpr([](double x) {
                            if ((x > (MIN_PRESSURE + PDD_DELTA)) &&
                                (x < (NORMAL_PRESSURE - PDD_DELTA)))
                              return 1.0;
                            return 0.0;
                          })
                          .array();

    auto pdd1_poly_vec = curves_info_->poly_coeffs()["PDD_POLY_VEC1"];
    auto pdd2_poly_vec = curves_info_->poly_coeffs()["PDD_POLY_VEC2"];

    auto vals =
        iso_junctions_.array() +
        connect_junctions_.array() *
            (case1_bool * (-PDD_SLOPE * demands_heads_vec_.array() *
                           variable_vec_->segment(0, nnodes_).array()) +
             case2_bool * (-demands_heads_vec_.array() *
                           (3 * pdd1_poly_vec[0] * pressure.array().pow(2) +
                            2 * pdd1_poly_vec[1] * pressure.array().pow(1) +
                            pdd1_poly_vec[2])) +
             case3_bool * (-demands_heads_vec_.array() *
                           (3 * pdd2_poly_vec[0] * pressure.array().pow(2) +
                            2 * pdd2_poly_vec[1] * pressure.array().pow(1) +
                            pdd2_poly_vec[2])) +
             case4_bool * (-PDD_SLOPE * demands_heads_vec_.array() *
                           variable_vec_->segment(0, nnodes_).array()) +
             case5_bool * (-0.5 * demands_heads_vec_.array() /
                           (NORMAL_PRESSURE - MIN_PRESSURE) *
                           ((pressure.array().abs() - MIN_PRESSURE) /
                            (NORMAL_PRESSURE - MIN_PRESSURE))
                               .abs()
                               .pow(-0.5)));
    //    std::cout << "====================================" << std::endl;
    //    std::cout << iso_junctions_.sum() << std::endl;
    //    std::cout <<case1_bool.sum()<<" "<< case2_bool.sum()<<" "<<
    //    case3_bool.sum()<<" "<< case4_bool.sum()<<" "<<case5_bool.sum()<<" "
    //    << std::endl; std::cout << vals.segment(0,nnodes_-source_idx_.size
    //    ()).sum() << std::endl;

    auto trip_d = sub_jac_trip_["jac_d"];
    for (int i = 0; i < nnodes_ - source_idx_.size(); ++i) {
      jac_->coeffRef(trip_d[i].row(), trip_d[i].col()) = vals[i];
    }
  }
}

// jac_f: Derivative of head-loss equation with respect to nodal head,
//              modified for powerpump
void pipenetwork::MatrixAssembler::update_jac_f() {
  auto trip_f = sub_jac_trip_["jac_f"];

  for (int i = 0; i < npumps_; ++i) {
    int pump_idx = npipes_ + i;
    auto link_flow = variable_vec_->array()[2 * nnodes_ + pump_idx];
    if (iso_links_[pump_idx] == 0) {
      auto pump = mesh_->links().at(pump_idx);
      auto pump_info = pump->link_info();
      if (pump_info["type"] == POWERPUMP) {
        jac_->coeffRef(trip_f[2 * pump_idx].row(), trip_f[2 * pump_idx].col()) =
            1000.0 * G * link_flow;
        jac_->coeffRef(trip_f[2 * pump_idx + 1].row(),
                       trip_f[2 * pump_idx + 1].col()) =
            -1000.0 * G * link_flow;
      }
    }
  }
}

// jac_g: Derivative of headloss equation with respect to pipe discharge,
//            for corresponding pipe, \frac{-1.852 \times 10.67 \times
//            length}{pow(pipe_roughness,1.852) \times pow(2radius,4.8704)}
//            \times pow(pipe_discharge,0.852).
void pipenetwork::MatrixAssembler::update_jac_g_pipe() {

  auto sign_array = (variable_vec_->segment(2 * nnodes_, npipes_))
                        .unaryExpr([](double x) {
                          if (x > 0) return 1.0;
                          return -1.0;
                        })
                        .array();  // get the sign of discharges

  // case 1, discharges that exceed the boundary of HW_Q2, use normal
  // hazen-william equation
  auto case1_bool = (variable_vec_->segment(2 * nnodes_, npipes_))
                        .unaryExpr([](double x) {
                          if (std::abs(x) > HW_Q2) return 1.0;
                          return 0.0;
                        })
                        .array();
  // case 2, discharges that fall in between HW_Q1 and HW_Q2, use polynomial
  // approximation
  auto case2_bool = (variable_vec_->segment(2 * nnodes_, npipes_))
                        .unaryExpr([](double x) {
                          if ((std::abs(x) < HW_Q2) && (std::abs(x) > HW_Q1))
                            return 1.0;
                          return 0.0;
                        })
                        .array();

  // case 3, discharges that are smaller than HW_Q1 , approximate 0 headloss in
  // this case
  auto case3_bool = (variable_vec_->segment(2 * nnodes_, npipes_))
                        .unaryExpr([](double x) {
                          if (std::abs(x) < HW_Q1) return 1.0;
                          return 0.0;
                        })
                        .array();
  auto discharge_abs_array =
      ((variable_vec_->segment(2 * nnodes_, npipes_)).array()).abs();
  auto head_diff_array = (headloss_mat_ * (variable_vec_->segment(0, nnodes_)))
                             .segment(0, npipes_)
                             .array();
  auto hw_poly_vec = curves_info_->poly_coeffs()["HW_POLY_VEC"];
  auto vals =
      (iso_links_.segment(0, npipes_)).array() +
      connect_links_.segment(0, npipes_).array() *
          (case1_bool * 1.852 *
               link_resistance_coeff_vec_.segment(0, npipes_).array() *
               discharge_abs_array.pow(.852) +
           case2_bool * link_resistance_coeff_vec_.segment(0, npipes_).array() *
               (3 * hw_poly_vec[0] * discharge_abs_array.pow(2) +
                2 * hw_poly_vec[1] * discharge_abs_array + 1 * hw_poly_vec[2]) +
           case3_bool * link_resistance_coeff_vec_.segment(0, npipes_).array() *
               HW_M);
  auto trip_g = sub_jac_trip_["jac_g"];
  for (int i = 0; i < npipes_; ++i) {
    jac_->coeffRef(trip_g[i].row(), trip_g[i].col()) = vals[i];
  }
}
// jac_g: Derivative of headloss equation with respect to pipe discharge,
//            for corresponding pump,
void pipenetwork::MatrixAssembler::update_jac_g_pump() {
  auto trip_g = sub_jac_trip_["jac_g"];
  // iterate through all the pumps
  for (int i = 0; i < npumps_; ++i) {
    int pump_idx = npipes_ + i;
    auto link_flow = (*variable_vec_)[2 * nnodes_ + pump_idx];
    if (iso_links_[pump_idx] == 1) {
      jac_->coeffRef(trip_g[pump_idx].row(), trip_g[pump_idx].col()) = 1;
    } else {
      auto pump = mesh_->links().at(pump_idx);
      auto pump_info = pump->link_info();
      auto nodes = pump->nodes();
      auto start_node_idx = node_id_map_.at(nodes.first->id());
      auto end_node_idx = node_id_map_.at(nodes.second->id());

      // head pump
      if (pump_info["type"] == HEADPUMP) {
        double pump_headgain;
        auto pump_curve = curves_info_->pump_curves().at(
            curves_info_->pump_int_str(pump_info.at("curve_name")));
        auto curve_coeff = pump_curve.head_curve_coefficients;

        if (curve_coeff[2] > 1) {
          auto line_coeff = pump_curve.line_param;
          if (link_flow >= line_coeff[0]) {
            jac_->coeffRef(trip_g[pump_idx].row(), trip_g[pump_idx].col()) =
                -curve_coeff[1] * curve_coeff[2] *
                std::pow(link_flow, (curve_coeff[2] - 1));

          }

          else {
            jac_->coeffRef(trip_g[pump_idx].row(), trip_g[pump_idx].col()) =
                PUMP_M;
          }
        } else {
          if (link_flow < PUMP_Q1) {
            jac_->coeffRef(trip_g[pump_idx].row(), trip_g[pump_idx].col()) =
                PUMP_M;
          } else if (link_flow < PUMP_Q2) {
            auto curve_poly_coeff = pump_curve.poly_coefficients;
            jac_->coeffRef(trip_g[pump_idx].row(), trip_g[pump_idx].col()) =
                3 * curve_poly_coeff[0] * std::pow(link_flow, 2) +
                2 * curve_poly_coeff[1] * link_flow + curve_poly_coeff[2];
          } else {
            jac_->coeffRef(trip_g[pump_idx].row(), trip_g[pump_idx].col()) =
                -curve_coeff[1] * curve_coeff[2] *
                std::pow(link_flow, (curve_coeff[2] - 1));
          }
        }
      }
      // power pump
      else if (pump_info["type"] == POWERPUMP) {

        auto head_diff_array =
            (headloss_mat_ * (variable_vec_->segment(0, nnodes_))).array();

        jac_->coeffRef(trip_g[pump_idx].row(), trip_g[pump_idx].col()) =
            (1000.0 * G * (*variable_vec_)[start_node_idx] -
             1000.0 * G * (*variable_vec_)[end_node_idx]);
      }
    }
  }
}

// jac_g: Derivative of headloss equation with respect to pipe discharge,
//            for corresponding valve,
void pipenetwork::MatrixAssembler::update_jac_g_valve() {
  auto head_diff_array =
      (headloss_mat_ * (variable_vec_->segment(0, nnodes_))).array();
  auto trip_g = sub_jac_trip_["jac_g"];
  // iterate through all the valves
  for (int i = 0; i < nvalves_; ++i) {
    int valve_idx = npipes_ + npumps_ + i;
    auto link_flow = (*variable_vec_)[2 * nnodes_ + valve_idx];

    if (iso_links_[valve_idx] == 1) {
      jac_->coeffRef(trip_g[valve_idx].row(), trip_g[valve_idx].col()) = 1;
    } else {

      auto valve = mesh_->links().at(valve_idx);
      auto valve_info = valve->link_info();
      auto nodes = valve->nodes();

      double val = 0;

      // PRV
      if (valve_info["type"] == PRVALVE) {
        if ((valve->link_status() == CLOSED)) {
          val = 1;
        } else if ((valve->link_status() == ACTIVE)) {
          val = 0;
        } else if ((valve->link_status() == OPEN)) {
          auto coeff = link_minor_loss_coeff_vec_[valve_idx];
          val = 2 * coeff * std::abs(link_flow);
        }
      }
      // FCV
      else if (valve_info["type"] == FCVALVE) {
        if ((valve->link_status() == CLOSED)) {
          val = 1;
        } else if ((valve->link_status() == ACTIVE)) {
          val = 1;
        } else if ((valve->link_status() == OPEN)) {
          auto coeff = link_minor_loss_coeff_vec_[valve_idx];
          val = 2 * coeff * std::abs(link_flow);
        }
      }
      // TCV
      else if (valve_info["type"] == TCVALVE) {
        if ((valve->link_status() == CLOSED)) {
          val = link_flow;
        } else if ((valve->link_status() == ACTIVE)) {
          auto coeff = link_resistance_coeff_vec_[valve_idx];
          val = 2 * coeff * std::abs(link_flow);
        } else if ((valve->link_status() == OPEN)) {
          auto coeff = link_minor_loss_coeff_vec_[valve_idx];
          val = 2 * coeff * std::abs(link_flow);
        }
      }
      jac_->coeffRef(trip_g[valve_idx].row(), trip_g[valve_idx].col()) = val;
    }
  }
}

// jac_h: Derivative of leak-pressure equation with respect to nodal head,
//              0 for inactive leaks, f(H-z) otherwise
void pipenetwork::MatrixAssembler::update_jac_h() {
  double m = 1e-11;
  Index leak_idx = 2 * nnodes_ + nlinks_;
  double val;
  auto trip_h = sub_jac_trip_["jac_h"];
  for (const auto& leak_id : leak_ids_) {
    auto i = leak_idx - 2 * nnodes_ - nlinks_;
    auto idx = node_id_map_[leak_id];
    auto p = (*variable_vec_)[idx] - elevations_[idx];

    // case 1, no pressure no leak
    if (p < m) {
      val = -m;
    }
    // case 2, around the boundary, use polynomial approximation
    else if (p < 1e-4) {
      auto leak_poly_coef = curves_info_->poly_coeffs()[leak_id];
      ;
      val = -3 * leak_poly_coef[0] * std::pow(p, 2) -
            2 * leak_poly_coef[1] * p - leak_poly_coef[2];
    }
    // case 3, normal leak equation
    else {
      val = -0.5 * LEAK_COEFF * leak_area_[leak_idx - 2 * nnodes_ - nlinks_] *
            std::sqrt(2 * G) * std::pow(p, -0.5);
    }
    jac_->coeffRef(trip_h[leak_idx - 2 * nnodes_ - nlinks_].row(),
                   trip_h[leak_idx - 2 * nnodes_ - nlinks_].col()) = val;
    ++leak_idx;
  }
}

void pipenetwork::MatrixAssembler::init_internal_graph() {

  std::map<std::pair<int, int>, int> n_links;
  std::vector<Eigen::Triplet<double>> graph_triplet;
  int val;
  // loop through links
  for (const auto& link : mesh_->links()) {
    auto nodes = link->nodes();
    auto start_node_idx = node_id_map_.at(nodes.first->id());
    auto end_node_idx = node_id_map_.at(nodes.second->id());
    // construct node link idx map
    if (node_link_id_map_.find(start_node_idx) == node_link_id_map_.end()) {
      std::vector<int> link_vec;
      node_link_id_map_[start_node_idx] = link_vec;
    }
    if (node_link_id_map_.find(end_node_idx) == node_link_id_map_.end()) {
      std::vector<int> link_vec;
      node_link_id_map_[end_node_idx] = link_vec;
    }
    node_link_id_map_[start_node_idx].emplace_back(link_id_map_[link->id()]);
    node_link_id_map_[end_node_idx].emplace_back(link_id_map_[link->id()]);

    // record number of links
    if (n_links.find(std::make_pair(start_node_idx, end_node_idx)) ==
        n_links.end()) {
      // not found
      n_links[std::make_pair(start_node_idx, end_node_idx)] = 0;
      n_links[std::make_pair(end_node_idx, start_node_idx)] = 0;
    }
    n_links[std::make_pair(start_node_idx, end_node_idx)] += 1;
    n_links[std::make_pair(end_node_idx, start_node_idx)] += 1;

    // put connectivity status as triplets
    val = 1;
    if (link->link_status() == CLOSED) {
      val = 0;
    }

    graph_triplet.emplace_back(start_node_idx, end_node_idx, val);
    graph_triplet.emplace_back(end_node_idx, start_node_idx, val);
  }
  internal_graph_.resize(nnodes_, nnodes_);
  internal_graph_.setFromTriplets(graph_triplet.begin(), graph_triplet.end());

  // create vector of number of connections
  nconnections_.resize(nnodes_);
  for (int i = 0; i < nnodes_; ++i) {
    nconnections_[i] = internal_graph_.outerIndexPtr()[i + 1] -
                       internal_graph_.outerIndexPtr()[i];
  }

  // init isolated junction and links arrays
  auto iso_nodes = get_isolated_nodes();
  auto iso_links = get_isolated_links(iso_nodes);
  iso_junctions_.setZero(nnodes_);
  iso_links_.setZero(nlinks_);
  for (auto& n : iso_nodes) {
    iso_junctions_[n] = 1;
  }
  for (auto& l : iso_links) {
    iso_links_[l] = 1;
  }
  Eigen::VectorXd ones_link;
  Eigen::VectorXd ones_nodes;
  ones_link.setOnes(nlinks_);
  ones_nodes.setOnes(nnodes_);
  connect_junctions_ = ones_nodes - iso_junctions_;
  connect_links_ = ones_link - iso_links_;
}

Eigen::VectorXd pipenetwork::MatrixAssembler::explore_nodes(int node_idx) {
  Eigen::VectorXd check_result(nnodes_);
  check_result.setZero(nnodes_);

  std::set<int> nodes_to_explore;
  check_result[node_idx] = 1;
  nodes_to_explore.emplace(node_idx);
  auto indptr = internal_graph_.outerIndexPtr();
  auto indices = internal_graph_.innerIndexPtr();
  auto data = internal_graph_.valuePtr();

  while (!nodes_to_explore.empty()) {
    auto node_being_explored = *nodes_to_explore.begin();
    nodes_to_explore.erase(nodes_to_explore.begin());
    //    std::cout<<"nodes being explore: "<<node_being_explored <<" "<<
    //    node_idx_map_[node_being_explored] <<std::endl;

    int nconnections = nconnections_[node_being_explored];
    int ndx = indptr[node_being_explored];
    // for all the connected nodes, set result to 1 and place them into the
    // searching queue
    for (int i = 0; i < nconnections; ++i) {
      if (data[ndx + i] != 0 && check_result[indices[ndx + i]] == 0) {
        check_result[indices[ndx + i]] = 1;
        nodes_to_explore.emplace(indices[ndx + i]);
      }
    }
  }
  return check_result;
}

std::vector<int> pipenetwork::MatrixAssembler::get_isolated_nodes() {
  Eigen::VectorXd check_result(nnodes_);
  check_result.setZero(nnodes_);
  std::vector<int> isolated_nodes;
  //  std::cout<<"souce number "<<source_idx_.size ()<<std::endl;
  for (const auto& source_idx : source_idx_) {
    check_result += explore_nodes(source_idx);
  }

  for (int i = 0; i < nnodes_; ++i) {
    if (check_result[i] == 0) {
      isolated_nodes.emplace_back(i);
      //      std::cout<< "iso node: "<< i << " " << node_idx_map_[i]
      //      <<std::endl;
    }
  }
  //    std::cout<<"iso junctions number "<<isolated_nodes.size ()<<std::endl;

  return isolated_nodes;
}

std::vector<int> pipenetwork::MatrixAssembler::get_isolated_links(
    const std::vector<int>& isolated_nodes) {
  std::vector<int> isolated_links;
  for (const auto& node : isolated_nodes) {
    for (const auto& link : node_link_id_map_[node]) {
      isolated_links.emplace_back(link);
    }
  }
  return isolated_links;
}