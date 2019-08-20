#include "hydralic_sim.h"
#include <iomanip>

bool pipenetwork::Hydralic_sim::run_simulation(double NR_tolerance,
                                               int max_nr_steps) {

  auto residual_vec = assembler_->residual_vector();
  auto variable_vec = assembler_->variable_vector();
  auto jac = assembler_->jac_matrix();

  solver_->assembled_matrices(jac, variable_vec, residual_vec);
  for (unsigned nr_iter = 0; nr_iter < max_nr_steps; ++nr_iter) {

    assembler_->assemble_residual();
    assembler_->update_jacobian();

    if (debug_) {
      // save the initial values into csv files for debugging
      if (nr_iter < 1) {
        std::ofstream outFile("../benchmarks/init_var_res.csv");
        std::ofstream outFile2("../benchmarks/init_jacob.csv");
        outFile << "variables"
                << ","
                << "residuals"
                << "\n";
        for (int i = 0; i < (*residual_vec).size(); ++i) {
          outFile << std::setprecision(12) << (*variable_vec).coeff(i) << ","
                  << (*residual_vec).coeff(i) << "\n";
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
      std::cout << "residual norm = " << residual_vec->norm() << std::endl;
      //                        std::cout << "Jac = " << std::endl
      //                                  << (*jac) << std::endl
      //                                  << std::endl
      //                                  << "residual = " << std::endl
      //                                  << (*residual_vec) << std::endl
      //                                  << std::endl
      //                                  << "variable = " << std::endl
      //                                  << (*variable_vec) << std::endl
      //                                  << std::endl;
    }

    bool issolved = solver_->solve();

    residual_norm_ = residual_vec->norm();
    if (residual_vec->norm() < NR_tolerance) {
      if (debug_) {
        std::ofstream outFile3("../benchmarks/final_var.csv");
        outFile3 << "variables"
                 << "\n";
        for (int i = 0; i < (*variable_vec).size(); ++i) {
          outFile3 << std::setprecision(12) << (*variable_vec).coeff(i) << "\n";
        }
        //        std::cout << "Final vairables " << (*variable_vec) <<
        //        std::endl;
      }
      std::string out_out_name = "../benchmarks/final_res";
      write_final_result(out_out_name, (*variable_vec));

      return true;
    }
  }
  return false;
}

pipenetwork::Hydralic_sim::Hydralic_sim(
    const std::string& filepath, const std::vector<double>& leak_diameters,
    bool pdd_mode, bool debug) {
  auto IO = std::make_shared<pipenetwork::Input>(filepath);

  // Mesh index
  const unsigned meshid = 9999;
  // Creat a mesh
  mesh_ = std::make_shared<pipenetwork::Mesh>(meshid);
    mesh_->create_mesh_from_inp(IO);
  // initialize discharges
    mesh_->iterate_over_links(std::bind(&pipenetwork::Link::update_sim_discharge,
                                  std::placeholders::_1,
                                  init_discharge_));  // initialze discharge
  // get curves information
  auto curves_info = IO->curve_info();
  assembler_ = std::make_shared<MatrixAssembler>(mesh_, curves_info, pdd_mode);
  solver_ = std::make_shared<Pardiso_unsym>(max_solver_steps_,
                                            inner_solver_tolerance_);
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
  auto nnodes = mesh_->nnodes ();
  auto nlinks = mesh_->nlinks ();

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
