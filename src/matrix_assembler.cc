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

  Index idx = 0, leak_idx = 0;
  for (const auto& node : mesh_->connect_nodes()) {
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

    ++idx;
  }

  // loop through links
  for (const auto& link : mesh_->links()) {
    Index index_pd = idx + nnodes_;
    // discharge vector
    variable_vec_->coeffRef(index_pd) = link->sim_discharge();
    // id map
    link_id_map_.emplace(link->id(), idx - nnodes_);

    auto info = link->link_info();
    double val = 0;
    switch (static_cast<pipenetwork::Link_type>(info["type"])) {
        // link resistence coefficients for pipes (Harzan-williams)
      case PIPE:
        val = HW_COEFF * std::pow(info["roughness"], -1.852) *
              std::pow(info["diameter"], -4.871) * info["length"];
        break;
      case TCVALVE:
        val = 8*info["setting"]/(G*std::pow(PI, 2)*std::pow(info["diameter"], 4));
        break;
    }
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
  // leak residual
  assemble_leak_residual();
}

// Residual for mass conservation equation (what comes in must go out)
void pipenetwork::MatrixAssembler::assemble_demand_head_residual() {
  if (!pdd_) {
    // demand for junctions
    residual_vec_->segment(nnodes_, nnodes_) =
        variable_vec_->segment(nnodes_, nnodes_) - demands_heads_vec_;

  } else {
    auto pressure = (variable_vec_->segment(0, nnodes_) - elevations_);

    // case 1, pressure smaller than min pressure, no water
    auto case1_bool = pressure
                          .unaryExpr([](double x) {
                            if (x < MIN_PRESSURE) return 1.0;
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
        case1_bool * (variable_vec_->segment(nnodes_, nnodes_).array()) +
        case2_bool *
            ((variable_vec_->segment(nnodes_, nnodes_).array()) -
             demands_heads_vec_.array() *
                 (pdd1_poly_vec[0] * pressure.array().pow(3) +
                  pdd1_poly_vec[1] * pressure.array().pow(2) +
                  pdd1_poly_vec[2] * pressure.array() + pdd1_poly_vec[3])) +
        case3_bool *
            ((variable_vec_->segment(nnodes_, nnodes_).array()) -
             demands_heads_vec_.array() *
                 (pdd2_poly_vec[0] * pressure.array().pow(3) +
                  pdd2_poly_vec[1] * pressure.array().pow(2) +
                  pdd2_poly_vec[2] * pressure.array() + pdd2_poly_vec[3])) +
        case4_bool * ((variable_vec_->segment(nnodes_, nnodes_).array()) -
                      demands_heads_vec_.array()) +
        case5_bool * ((variable_vec_->segment(nnodes_, nnodes_).array()) -
                      demands_heads_vec_.array() *
                          ((pressure.array().abs() - MIN_PRESSURE) /
                           (NORMAL_PRESSURE - MIN_PRESSURE))
                              .pow(0.5));
    //    std::cout<<pressure.array() <<std::endl;
  }
  // correct residuals for sources (head for reservoir/tanks)
  for (const auto& idx : source_idx_) {
    (*residual_vec_)[idx + nnodes_] =
        (*variable_vec_)[idx] - demands_heads_vec_[idx];
  }
}

// Residual for leak from holes. piece-wise function, use for loop... (number of
// leak nodes is usually small too)
void pipenetwork::MatrixAssembler::assemble_leak_residual() {
  double m = 1e-11;
  Index leak_idx = 2 * nnodes_ + nlinks_;
  for (const auto& leak_id : leak_ids_) {
    auto i = leak_idx - 2 * nnodes_ - nlinks_;
    auto idx = node_id_map_[leak_id];
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
      case1_bool *
          (sign_array * link_resistance_coeff_vec_.segment(0, npipes_).array() *
               discharge_abs_array.pow(1.852) -
           head_diff_array)

      +
      case2_bool *
          (sign_array * link_resistance_coeff_vec_.segment(0, npipes_).array() *
               (hw_poly_vec[0] * discharge_abs_array.pow(3) +
                hw_poly_vec[1] * discharge_abs_array.pow(2) +
                hw_poly_vec[2] * discharge_abs_array + hw_poly_vec[3]) -
           head_diff_array) +
      case3_bool *
          (sign_array * link_resistance_coeff_vec_.segment(0, npipes_).array() *
               HW_M * discharge_abs_array -
           head_diff_array);
}

// TODO: test pump part
void pipenetwork::MatrixAssembler::assemble_headloss_residual_pump() {
  // iterate through all the pumps
  for (int i = 0; i < npumps_; ++i) {
    int pump_idx = npipes_ + i;
    auto link_flow = (*variable_vec_)[2 * nnodes_ + pump_idx];

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
          pump_headgain = curve_coeff[0] -
                          curve_coeff[1] * std::pow(link_flow, curve_coeff[2]);

        }

        else {
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
          pump_headgain = curve_coeff[0] -
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
  // Jac_d: pressure-demand equation
  update_jac_d();
  // Jac_f: for power pump only
  update_jac_f();
  // jac_g: headloss equation (harzen-william)
  update_jac_g_pipe();
  update_jac_g_pump();
  // jac_h: leak equation
  update_jac_h();
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
                            if (x < MIN_PRESSURE) return 1.0;
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
        case1_bool *
            (-PDD_SLOPE * variable_vec_->segment(nnodes_, nnodes_).array()) +
        case2_bool * (-demands_heads_vec_.array() *
                      (3 * pdd1_poly_vec[0] * pressure.array().pow(2) +
                       2 * pdd1_poly_vec[1] * pressure.array().pow(1) +
                       pdd1_poly_vec[2])) +
        case3_bool * (-demands_heads_vec_.array() *
                      (3 * pdd2_poly_vec[0] * pressure.array().pow(2) +
                       2 * pdd2_poly_vec[1] * pressure.array().pow(1) +
                       pdd2_poly_vec[2])) +
        case4_bool *
            (-PDD_SLOPE * variable_vec_->segment(nnodes_, nnodes_).array()) +
        case5_bool * (-0.5 * demands_heads_vec_.array() *
                      ((pressure.array().abs() - MIN_PRESSURE) /
                       (NORMAL_PRESSURE - MIN_PRESSURE))
                          .pow(-0.5));

    auto trip_d = sub_jac_trip_["jac_d"];
    for (int i = 0; i < nnodes_; ++i) {
      if (trip_d[i].value() == 0)
        jac_->coeffRef(trip_d[i].row(), trip_d[i].col()) = vals[i];
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
      case1_bool * 1.852 *
          link_resistance_coeff_vec_.segment(0, npipes_).array() *
          discharge_abs_array.pow(.852) +
      case2_bool * link_resistance_coeff_vec_.segment(0, npipes_).array() *
          (3 * hw_poly_vec[0] * discharge_abs_array.pow(2) +
           2 * hw_poly_vec[1] * discharge_abs_array + 1 * hw_poly_vec[2]) +
      case3_bool * link_resistance_coeff_vec_.segment(0, npipes_).array() *
          HW_M;
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
// jac_f: Derivative of head-loss equation with respect to nodal head,
//              modified for powerpump
void pipenetwork::MatrixAssembler::update_jac_f() {
  auto trip_f = sub_jac_trip_["jac_f"];

  for (int i = 0; i < npumps_; ++i) {
    int pump_idx = npipes_ + i;
    auto link_flow = variable_vec_->array()[2 * nnodes_ + pump_idx];

    auto pump = mesh_->links().at(pump_idx);
    auto pump_info = pump->link_info();
    if (pump_info["type"] == POWERPUMP) {
      jac_->coeffRef(trip_f[2 * pump_idx].row(), trip_f[2 * pump_idx].col()) =
          1000.0 * G * link_flow;
      jac_->coeffRef(trip_f[2 * pump_idx + 1].row(),
                     trip_f[2 * pump_idx + 1].col()) = -1000.0 * G * link_flow;
    }
  }
}
