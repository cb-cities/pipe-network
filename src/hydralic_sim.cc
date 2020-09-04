#include "hydralic_sim.h"

bool pipenetwork::Hydralic_sim::run_simulation(double NR_tolerance,
                                               int max_nr_steps) {

  solver_->assembled_matrices(assembler_);
  for (unsigned nr_iter = 0; nr_iter < max_nr_steps; ++nr_iter) {

    assembler_->system_update();
    init_variable_ = assembler_->variable_vector();

    auto x_diff = solver_->solve();
    line_search_update(x_diff);
    auto res_norm = assembler_->residual_vector().norm();

    if (debug_) {
      std::cout << "niter = " << nr_iter << std::endl;
      std::cout << "residual norm = " << assembler_->residual_vector().norm()
                << std::endl;
      //      std::cout << "variable " << assembler_->variable_vector() <<
      //      std::endl; std::cout << "residual " <<
      //      assembler_->residual_vector() << std::endl;
    }

    if (res_norm < NR_tolerance) {
      return true;
    }
  }
  return false;
}

void pipenetwork::Hydralic_sim::line_search_update(
    const Eigen::VectorXd& x_diff) {
  double alpha = 1.0;  // start line search with the original x_diff

  auto old_res = assembler_->residual_vector().norm();
  Eigen::VectorXd& variables = assembler_->variable_vector();

  for (int bt_iter = 0; bt_iter < bt_max_iter_; ++bt_iter) {

    variables = init_variable_.array() - alpha * x_diff.array();
    assembler_->assemble_residual();  // update residual from the new variable
                                      // vector
    auto new_res = assembler_->residual_vector().norm();
    if (new_res < (1.0 - 0.0001 * alpha) * old_res) {
      break;  // line search finished
    } else {
      alpha = alpha * bt_roh_;  // continue line search
    }
  }
}

pipenetwork::Hydralic_sim::Hydralic_sim(
    const std::shared_ptr<pipenetwork::Mesh>& mesh,
    std::shared_ptr<pipenetwork::Curves>& curves_info, bool pdd_mode,
    bool debug) {
  mesh_ = mesh;
  assembler_ = std::make_shared<linear_system::MatrixAssembler>(
      mesh, curves_info, pdd_mode);
  // choose solver type based on network size
  Index nnodes = mesh_->nodes()->nnodes();
  if (nnodes < nthre_) {
    std::shared_ptr<linear_system::Solver> solve_ptr(
        Factory<linear_system::Solver>::instance()->create("LU"));
    std::cout << "small network, using the LU Solver " << std::endl;
    solver_ = solve_ptr;
  } else {
    std::shared_ptr<linear_system::Solver> solve_ptr(
        Factory<linear_system::Solver>::instance()->create("mkl_pardiso"));
    std::cout << "large network, using the mkl_pardiso Solver " << std::endl;
    solver_ = solve_ptr;
  }

  debug_ = debug;
}

void pipenetwork::Hydralic_sim::update_mesh() {
  auto nodes = mesh_->nodes();
  auto links = mesh_->links();
  auto leak_nids = mesh_->leak_nids();

  auto variables = assembler_->variable_vector();
  auto nnodes = nodes->nnodes();
  auto nlinks = links->nlinks();

  // junctions
  auto junction_map = nodes->junctions();
  Index leak_count = 0;
  for (auto& index_junc : junction_map) {
    auto nid = index_junc.first;
    auto junction = index_junc.second;
    auto sim_head = variables[nid];
    auto sim_demand = variables[nnodes + nid];
    junction->head = sim_head;
    junction->demand = sim_demand;
    bool leak_node =
        std::find(leak_nids.begin(), leak_nids.end(), nid) != leak_nids.end();
    if (leak_node) {
      junction->leak_discharge = variables[2 * nnodes + nlinks + leak_count];
      leak_count++;
    }
  }

  // srcs
  auto res_map = nodes->reservoirs();
  for (auto& index_res : res_map) {
    auto nid = index_res.first;
    auto reservoir = index_res.second;
    auto sim_discharge = variables[nnodes + nid];
    reservoir->discharge = sim_discharge;
  }

  // pipes
  auto pipe_map = links->pipes();
  for (auto& index_pipe : pipe_map) {
    auto lid = index_pipe.first;
    auto pipe = index_pipe.second;
    auto sim_flowrate = variables[2 * nnodes + lid];
    pipe->flowrate = sim_flowrate;
  }
}
