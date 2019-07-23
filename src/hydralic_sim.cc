#include "hydralic_sim.h"

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
          outFile << (*variable_vec).coeff(i) << "," << (*residual_vec).coeff(i)
                  << "\n";
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
            outFile2 << it.row() << "," << it.col() << "," << it.value()
                     << "\n";
          }
      }
      std::cout << "niter = " << nr_iter << std::endl;
      std::cout << "residual norm = " << residual_vec->norm() << std::endl;
      //                  std::cout << "Jac = " << std::endl
      //                            << (*jac) << std::endl
      //                            << std::endl
      //                            << "residual = " << std::endl
      //                            << (*residual_vec) << std::endl
      //                            << std::endl
      //                            << "variable = " << std::endl
      //                            << (*variable_vec) << std::endl
      //                            << std::endl;
    }

    bool issolved = solver_->solve();

    residual_norm_ = residual_vec->norm();
        if (residual_vec->norm() < NR_tolerance) {
          if (debug_) {
            std::ofstream outFile3("../benchmarks/final_var.csv");
            outFile3 << "variables"
                     << "\n";
            for (int i = 0; i < (*variable_vec).size(); ++i) {
              outFile3 << (*variable_vec).coeff(i) << "\n";
            }
            std::cout << "Final vairables " << (*variable_vec) << std::endl;
          }

          return true;
        }
  }
  return false;
}

pipenetwork::Hydralic_sim::Hydralic_sim(
    const std::string& filepath, const std::vector<double>& leak_diameters,
    bool pdd_mode, bool debug) {
  auto IO = std::make_unique<pipenetwork::Input>(filepath);

  // Mesh index
  const unsigned meshid = 9999;
  // Creat a mesh
  auto m = std::make_shared<pipenetwork::Mesh>(meshid);
  m->create_junctions(IO->junction_properties());
  m->create_reservoirs(IO->reservoir_properties());
  std::vector<pipenetwork::Pipe_prop> pipe_props = IO->pipe_properties();
  m->create_pipes(pipe_props);
  // initialize discharges
  m->iterate_over_links(std::bind(&pipenetwork::Link::update_sim_discharge,
                                  std::placeholders::_1,
                                  init_discharge_));  // initialze discharge
  assembler_ = std::make_shared<MatrixAssembler>(m, pdd_mode);
  solver_ = std::make_shared<Pardiso_unsym>(max_solver_steps_,
                                            inner_solver_tolerance_);
  debug_ = debug;
  // print mesh summary if on debug mode
  if (debug_) m->print_summary();
}
