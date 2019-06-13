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
    //
    //    std::cout << "niter = " << nr_iter << std::endl;
    //    std::cout << "residual norm = " << residual_vec->norm() << std::endl;
    //    std::cout << "Jac = " << std::endl
    //              << (*jac) << std::endl
    //              << std::endl
    //              << "residual = " << std::endl
    //              << (*residual_vec) << std::endl
    //              << std::endl
    //              << "variable = " << std::endl
    //              << (*variable_vec) << std::endl
    //              << std::endl;

    bool issolved = solver_->solve();
    if (residual_vec->norm() < NR_tolerance) {
      residual_norm_ = residual_vec->norm();
      return true;
    }
  }
  residual_norm_ = residual_vec->norm();
  return false;
}

pipenetwork::Hydralic_sim::Hydralic_sim(
    const std::string & filepath, const std::vector<double>& leak_diameters,
    bool pdd_mode) {
  auto IO = std::make_unique<pipenetwork::Input>(filepath);

  // Mesh index
  const unsigned meshid = 9999;
  // Creat a mesh
  auto m = std::make_shared<pipenetwork::Mesh>(meshid);
  m->create_junctions(IO->junction_ids(), IO->junction_elevations(),
                      IO->junction_demands(), leak_diameters);
  m->create_reservoirs(IO->reservoir_ids(), IO->reservoir_heads());
  m->create_pipes(IO->pipe_ids(), IO->pipe_nodes_ids(), IO->pipe_length(),
                  IO->pipe_diameters(), IO->pipe_roughness(),
                  IO->pipe_status());
  //  m->print_summary();
  // initialize discharges
  m->iterate_over_links(std::bind(&pipenetwork::Link::update_sim_discharge,
                                  std::placeholders::_1,
                                  init_discharge_));  // initialze discharge
  assembler_ = std::make_shared<MatrixAssembler>(m, pdd_mode);
  solver_ =
      std::make_shared<EigenGMRES>(max_solver_steps_, inner_solver_tolerance_);
}
