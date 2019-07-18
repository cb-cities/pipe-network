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
          leak_ids_.emplace_back(node.second->id());
          leak_area_.emplace_back(info["leak_area"]);
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
    variable_vec_->coeffRef(index_pd) = link.second->sim_discharge();
    // id map
    link_id_map_.emplace(link.second->id(), idx - nnodes_);
    // link resistence coefficients for pipes
    auto info = link.second->link_info();
    switch (static_cast<int>(info["type"])) {
      case PIPE:
        link_resistance_coeff_vec_[idx - nnodes_] =
            HW_COEFF * std::pow(info["roughness"], -1.852) *
            std::pow(info["diameter"], -4.871) * info["length"];
        break;
    }
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
    auto out_node_id = link.second->nodes().first->id();
    auto in_node_id = link.second->nodes().second->id();

    update_balance.emplace_back(node_id_map_[out_node_id],
                                link_id_map_[link.second->id()], -1);
    update_balance.emplace_back(node_id_map_[in_node_id],
                                link_id_map_[link.second->id()], 1);

    update_jacb.emplace_back(node_id_map_[out_node_id],
                             link_id_map_[link.second->id()] + 2 * nnodes_, -1);
    update_jacb.emplace_back(node_id_map_[in_node_id],
                             link_id_map_[link.second->id()] + 2 * nnodes_, 1);

    update_headloss.emplace_back(link_id_map_[link.second->id()],
                                 node_id_map_[out_node_id], 1);
    update_headloss.emplace_back(link_id_map_[link.second->id()],
                                 node_id_map_[in_node_id], -1);

    update_jacf.emplace_back(link_id_map_[link.second->id()] + 2 * nnodes_,
                             node_id_map_[out_node_id], -1);
    update_jacf.emplace_back(link_id_map_[link.second->id()] + 2 * nnodes_,
                             node_id_map_[in_node_id], 1);
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
  assemble_headloss_residual();
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
    residual_vec_->segment(nnodes_, nnodes_) =
        case1_bool * (variable_vec_->segment(nnodes_, nnodes_).array()) +
        case2_bool *
            ((variable_vec_->segment(nnodes_, nnodes_).array()) -
             demands_heads_vec_.array() *
                 (PDD_POLY_VEC1[0] * pressure.array().pow(3) +
                  PDD_POLY_VEC1[1] * pressure.array().pow(2) +
                  PDD_POLY_VEC1[2] * pressure.array() + HW_POLY_VEC[3])) +
        case3_bool *
            ((variable_vec_->segment(nnodes_, nnodes_).array()) -
             demands_heads_vec_.array() *
                 (PDD_POLY_VEC2[0] * pressure.array().pow(3) +
                  PDD_POLY_VEC2[1] * pressure.array().pow(2) +
                  PDD_POLY_VEC2[2] * pressure.array() + HW_POLY_VEC[3])) +
        case4_bool * ((variable_vec_->segment(nnodes_, nnodes_).array()) -
                      demands_heads_vec_.array()) +
        case5_bool *
            ((variable_vec_->segment(nnodes_, nnodes_).array()) -
             demands_heads_vec_.array() * ((pressure.array() - MIN_PRESSURE) /
                                           (NORMAL_PRESSURE - MIN_PRESSURE))
                                              .pow(0.5));

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
      auto leak_poly_coef = compute_leak_poly_coef(leak_area_[i]);
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
void pipenetwork::MatrixAssembler::assemble_headloss_residual() {
  auto sign_array = (variable_vec_->segment(2 * nnodes_, nlinks_))
                        .unaryExpr([](double x) {
//                            return (x < 0) ? -1 : 1;
                            if (x > 0) return 1.0;
                            return -1.0;
                        })
                        .array();  // get the sign of discharges

  // case 1, discharges that exceed the boundary of HW_Q2, use normal
  // hazen-william equation
  auto case1_bool = (variable_vec_->segment(2 * nnodes_, nlinks_))
                        .unaryExpr([](double x) {
                            if (std::abs(x) > HW_Q2) return 1.0;
                            return 0.0;


                        })
                        .array();
  // case 2, discharges that fall in between HW_Q1 and HW_Q2, use polynomial
  // approximation
  auto case2_bool = (variable_vec_->segment(2 * nnodes_, nlinks_))
                        .unaryExpr([](double x) {
                          if ((std::abs(x) < HW_Q2) && (std::abs(x) > HW_Q1))
                            return 1.0;
                          return 0.0;
                        })
                        .array();
  // case 3, discharges that are smaller than HW_Q1 , approximate 0 headloss in
  // this case
  auto case3_bool = (variable_vec_->segment(2 * nnodes_, nlinks_))
                        .unaryExpr([](double x) {
                          if (std::abs(x) < HW_Q1) return 1.0;
                          return 0.0;
                        })
                        .array();
  auto discharge_abs_array =
      ((variable_vec_->segment(2 * nnodes_, nlinks_)).array()).abs();
  auto head_diff_array =
      (headloss_mat_ * (variable_vec_->segment(0, nnodes_))).array();

  residual_vec_->segment(2 * nnodes_, nlinks_) =
      case1_bool * (sign_array * link_resistance_coeff_vec_.array() *
                        discharge_abs_array.pow(1.852) -
                    head_diff_array)

      + case2_bool *
            (sign_array * link_resistance_coeff_vec_.array() *
                 (HW_POLY_VEC[0] * discharge_abs_array.pow(3) +
                  HW_POLY_VEC[1] * discharge_abs_array.pow(2) +
                  HW_POLY_VEC[2] * discharge_abs_array + HW_POLY_VEC[3]) -
             head_diff_array) +
      case3_bool * (sign_array * link_resistance_coeff_vec_.array() * HW_M *
                        discharge_abs_array -
                    head_diff_array);
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
  // jac_g: headloss equation (harzen-william)
  update_jac_g();
  // jac_h: leak equation
  update_jac_h();
}

// jac_d: Derivative of demand_pressure equation with respect to nodal head,
//            0 for junctions, 1 for reservoir/tank of demand driven model,
//            nonlinear for pressure_demand model.
// TODO: pressure demand part
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

    auto vals =
        case1_bool *
            (-PDD_SLOPE * variable_vec_->segment(nnodes_, nnodes_).array()) +
        case2_bool * (-demands_heads_vec_.array() *
                      (3 * PDD_POLY_VEC1[0] * pressure.array().pow(2) +
                       2 * PDD_POLY_VEC1[1] * pressure.array().pow(1) +
                       PDD_POLY_VEC1[2])) +
        case3_bool * (-demands_heads_vec_.array() *
                      (3 * PDD_POLY_VEC2[0] * pressure.array().pow(2) +
                       2 * PDD_POLY_VEC2[1] * pressure.array().pow(1) +
                       PDD_POLY_VEC2[2])) +
        case4_bool *
            (-PDD_SLOPE * variable_vec_->segment(nnodes_, nnodes_).array()) +
        case5_bool * (-0.5 * demands_heads_vec_.array() *
                      ((pressure.array() - MIN_PRESSURE) /
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
void pipenetwork::MatrixAssembler::update_jac_g() {

  auto sign_array = (variable_vec_->segment(2 * nnodes_, nlinks_))
                        .unaryExpr([](double x) {
                            if (x > 0) return 1.0;
                            return -1.0;
                        })
                        .array();  // get the sign of discharges

  // case 1, discharges that exceed the boundary of HW_Q2, use normal
  // hazen-william equation
  auto case1_bool = (variable_vec_->segment(2 * nnodes_, nlinks_))
                        .unaryExpr([](double x) {
                          if (std::abs(x) > HW_Q2) return 1.0;
                          return 0.0;
                        })
                        .array();
  // case 2, discharges that fall in between HW_Q1 and HW_Q2, use polynomial
  // approximation
  auto case2_bool = (variable_vec_->segment(2 * nnodes_, nlinks_))
                        .unaryExpr([](double x) {
                          if ((std::abs(x) < HW_Q2) && (std::abs(x) > HW_Q1))
                            return 1.0;
                          return 0.0;
                        })
                        .array();

  // case 3, discharges that are smaller than HW_Q1 , approximate 0 headloss in
  // this case
  auto case3_bool = (variable_vec_->segment(2 * nnodes_, nlinks_))
                        .unaryExpr([](double x) {
                          if (std::abs(x) < HW_Q1) return 1.0;
                          return 0.0;
                        })
                        .array();
  auto discharge_abs_array =
      ((variable_vec_->segment(2 * nnodes_, nlinks_)).array()).abs();
  auto head_diff_array =
      (headloss_mat_ * (variable_vec_->segment(0, nnodes_))).array();

  auto vals =
      case1_bool * 1.852 * link_resistance_coeff_vec_.array() *
          discharge_abs_array.pow(.852) +
      case2_bool * link_resistance_coeff_vec_.array() *
          (3 * HW_POLY_VEC[0] * discharge_abs_array.pow(2) +
           2 * HW_POLY_VEC[1] * discharge_abs_array + 1 * HW_POLY_VEC[2]) +
      case3_bool * link_resistance_coeff_vec_.array() * HW_M;

  auto trip_g = sub_jac_trip_["jac_g"];
  for (int i = 0; i < nlinks_; ++i) {
    jac_->coeffRef(trip_g[i].row(), trip_g[i].col()) = vals[i];
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
      auto leak_poly_coef = compute_leak_poly_coef(leak_area_[i]);
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

Eigen::VectorXd pipenetwork::MatrixAssembler::compute_poly_coefficients(
    const std::array<double, 2>& x, const std::array<double, 2>& f,
    const std::array<double, 2>& df) const {
  Eigen::VectorXd ret(4);
  double a =
      (2 * (f[0] - f[1]) - (x[0] - x[1]) * (df[1] + df[0])) /
      (std::pow(x[1], 3) - std::pow(x[0], 3) + 3 * x[0] * x[1] * (x[0] - x[1]));
  double b =
      (df[0] - df[1] + 3 * ((std::pow(x[1], 2) - std::pow(x[0], 2)) * a)) /
      (2 * (x[0] - x[1]));
  double c = df[1] - 3 * std::pow(x[1], 2) * a - 2 * x[1] * b;
  double d = f[1] - std::pow(x[1], 3) * a - std::pow(x[1], 2) * b - x[1] * c;
  ret << a, b, c, d;
  return ret;
}

Eigen::VectorXd pipenetwork::MatrixAssembler::compute_leak_poly_coef(
    double leak_area) const {
  double x1 = 0.0;
  double f1 = 0.0;
  double df1 = 1.0e-11;
  double x2 = 1e-4;
  double f2 = LEAK_COEFF * leak_area * std::pow((2 * G * x2), 0.5);
  double df2 = 0.5 * LEAK_COEFF * leak_area * std::pow((2 * G), 0.5) *
               std::pow(x2, -0.5);

  return compute_poly_coefficients({x1, x2}, {f1, f2}, {df1, df2});
}
