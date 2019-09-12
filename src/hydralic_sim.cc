#include "hydralic_sim.h"
#include <iomanip>

bool pipenetwork::Hydralic_sim::run_simulation(double NR_tolerance,
                                               int max_nr_steps,
                                               bool line_search,
                                               std::string output_path) {

  residuals_ = assembler_->residual_vector();
  variables_ = assembler_->variable_vector();
  auto jac = assembler_->jac_matrix();

  solver_->assembled_matrices(jac, variables_, residuals_);
  for (unsigned nr_iter = 0; nr_iter < max_nr_steps; ++nr_iter) {

    assembler_->assemble_residual();
    assembler_->update_jacobian();
    original_variable_ = *variables_;

    if (debug_) {
      // save the initial values into csv files for debugging
      if (nr_iter < 1) {
        std::ofstream outFile("../benchmarks/init_var_res.csv");
        std::ofstream outFile2("../benchmarks/init_jacob.csv");
        outFile << "variables"
                << ","
                << "residuals"
                << "\n";
        for (int i = 0; i < (*residuals_).size(); ++i) {
          outFile << std::setprecision(12) << (*variables_).coeff(i) << ","
                  << (*residuals_).coeff(i) << "\n";
        }
        outFile2 << "row"
                 << ","
                 << "col"
                 << ","
                 << "val"
                 << "\n";
        for (int k = 0; k < (*jac).outerSize(); ++k)
          for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(
                   (*jac), k);
               it; ++it) {
            outFile2 << std::setprecision(12) << it.row() << "," << it.col()
                     << "," << it.value() << "\n";
          }
      }
      std::cout << "niter = " << nr_iter << std::endl;
      std::cout << "residual norm = " << residuals_->norm() << std::endl;
      //                                    std::cout << "Jac = " << std::endl
      //                                              << (*jac) << std::endl
      //                                              << std::endl
      //                                              << "residual = " <<
      //                                              std::endl
      //                                              << (*residual_vec) <<
      //                                              std::endl
      //                                              << std::endl
      //                                              << "variable = " <<
      //                                              std::endl
      //                                              << (*variable_vec) <<
      //                                              std::endl
      //                                              << std::endl;
    }
    //    std::cout << "start solving" << std::endl;
    auto x_diff = solver_->solve();

    //    std::cout << "start line search" << std::endl;
    if (line_search) {
      line_search_func(x_diff);
    } else {
      (*variables_) = original_variable_.array() - x_diff.array();
    }
    residual_norm_ = residuals_->norm();
    if (residuals_->norm() < NR_tolerance) {
      auto path_name = output_path + mesh_->id();
      write_final_result(path_name, (*variables_));
      return true;
    }
  }
  return false;
}
void pipenetwork::Hydralic_sim::line_search_func(
    const Eigen::VectorXd& x_diff) {
  double alpha = 1.0;
  auto old_res = residuals_->norm();
  for (int bt_iter = 0; bt_iter < bt_max_iter_; ++bt_iter) {
    (*variables_) = original_variable_.array() - alpha * x_diff.array();
    assembler_->assemble_residual();
    auto new_res = residuals_->norm();
    if (new_res < (1.0 - 0.0001 * alpha) * old_res) {
      break;
    } else {
      alpha = alpha * bt_roh_;
    }
  }
  //  throw std::runtime_error("Line search failed!");
}

pipenetwork::Hydralic_sim::Hydralic_sim(const std::string& filepath,
                                        const std::string& mesh_name,
                                        bool pdd_mode,
                                        const std::string& solver_name,
                                        bool debug) {
  auto IO = std::make_shared<pipenetwork::Input>(filepath);
  // Creat a mesh
  std::string mesh_id = mesh_name;
  mesh_ = std::make_shared<pipenetwork::Mesh>(mesh_id);
  mesh_->create_mesh_from_inp(IO);
  // initialize discharges
  mesh_->iterate_over_links(std::bind(&pipenetwork::Link::update_sim_discharge,
                                      std::placeholders::_1,
                                      init_discharge_));  // initialze discharge
  // get curves information
  auto curves_info = IO->curve_info();
  assembler_ = std::make_shared<MatrixAssembler>(mesh_, curves_info, pdd_mode);

  //  solver_ = std::make_shared<Mkl_unsym>();
  std::shared_ptr<Solver> solve_ptr(
      Factory<Solver>::instance()->create(solver_name));
  solver_ = solve_ptr;
  debug_ = debug;
  // print mesh summary if on debug mode
  if (debug_) mesh_->print_summary();
}

pipenetwork::Hydralic_sim::Hydralic_sim(int syn_size, bool pdd_mode,
                                        const std::string& solver_name,
                                        bool debug) {
  auto IO = std::make_shared<pipenetwork::Input>(syn_size);

  // Mesh index
  std::string meshid = "syn_net_" + std::to_string(syn_size);
  // Creat a mesh
  mesh_ = std::make_shared<pipenetwork::Mesh>(meshid);
  mesh_->create_mesh_from_inp(IO);
  auto Output = std::make_shared<pipenetwork::Output>(
      mesh_, "../benchmarks/write_test_large.inp");
  // initialize discharges
  mesh_->iterate_over_links(std::bind(&pipenetwork::Link::update_sim_discharge,
                                      std::placeholders::_1,
                                      init_discharge_));  // initialze discharge
  // get curves information
  auto curves_info = IO->curve_info();
  assembler_ = std::make_shared<MatrixAssembler>(mesh_, curves_info, pdd_mode);
  std::shared_ptr<Solver> solve_ptr(
      Factory<Solver>::instance()->create(solver_name));
  solver_ = solve_ptr;
  debug_ = debug;
  // print mesh summary if on debug mode
  if (debug_) mesh_->print_summary();
}

void pipenetwork::Hydralic_sim::write_final_result(
    const std::string& output_path, const Eigen::VectorXd& var) {
  std::ofstream outnode(output_path + "_nodes.csv");
  std::ofstream outlink(output_path + "_links.csv");
  outnode << "node_id"
          << ","
          << "head"
          << ","
          << "demand"
          << "\n";
  outlink << "link_id"
          << ","
          << "flowrate"
          << "\n";
  auto node_map = assembler_->node_idx_map();
  auto link_map = assembler_->link_idx_map();
  auto nnodes = mesh_->nnodes();
  auto nlinks = mesh_->nlinks();

  std::string node_id, link_id;
  int demand_idx;
  double head, demand, flow_rate;
  for (int i = 0; i < nnodes; ++i) {
    demand_idx = nnodes + i;
    node_id = node_map.at(i);
    head = var[i];
    demand = var[demand_idx];
    outnode << std::setprecision(12) << node_id << "," << head << "," << demand
            << "\n";
  }
  for (int i = 0; i < nlinks; ++i) {
    link_id = link_map.at(i);
    flow_rate = var[2 * nnodes + i];
    outlink << std::setprecision(12) << link_id << "," << flow_rate << "\n";
  }
}
