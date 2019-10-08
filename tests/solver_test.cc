#include "catch.hpp"

#include "curves.h"
#include "factory.h"
#include "input.h"
#include "matrix_assembler.h"
#include "solver.h"

// Check solver class
TEST_CASE("Solver is checked", "[Solver]") {

  // Tolerance
  const double tolerance = 1.e-6;

  // Mesh index
  std::string meshid = "Matrix test mesh";

  // Creat a mesh
  auto mesh = std::make_shared<pipenetwork::Mesh>(meshid);

  // Create a curves info object
  auto curves_info = std::make_shared<pipenetwork::Curves>();

  std::vector<std::string> junction_ids{"10", "11", "12"};
  std::vector<double> elevations{216.408, 216.408, 213.36};
  std::vector<double> demands{0, 9.464e-03, 9.464e-03};
  std::vector<double> leak_diameters{0, 0.1, 0};

  std::vector<pipenetwork::Junction_prop> junc_props;
  for (int i = 0; i < elevations.size(); ++i) {
    pipenetwork::Junction_prop junc_prop;
    junc_prop.id = junction_ids[i];
    junc_prop.elevation = elevations[i];
    junc_prop.demand = demands[i];
    junc_prop.leak_diameter = leak_diameters[i];

    junc_props.emplace_back(junc_prop);
  }

  mesh->create_junctions(junc_props);

  std::vector<std::string> res_ids{"13"};
  std::vector<double> heads{3.048e+02};

  std::vector<pipenetwork::Reservoir_prop> res_props;
  for (int i = 0; i < res_ids.size(); ++i) {
    pipenetwork::Reservoir_prop res_prop;
    res_prop.id = res_ids[i];
    res_prop.head = heads[i];

    res_props.emplace_back(res_prop);
  }

  mesh->create_reservoirs(res_props);

  std::vector<std::string> pipe_ids{"10", "11", "12", "13"};
  std::vector<std::pair<std::string, std::string>> nodeids{
      std::make_pair("13", "10"), std::make_pair("10", "11"),
      std::make_pair("11", "12"), std::make_pair("10", "12")};
  const std::vector<double> length{3209.5440000000003, 3209.5440000000003,
                                   1609.344, 1609.344};
  const std::vector<double> diameter{0.5588, 0.4572, 0.35559999999999997,
                                     0.254};
  const std::vector<double> roughness{100, 100, 100, 100};
  const std::vector<pipenetwork::Link_status> status{
      pipenetwork::OPEN, pipenetwork::OPEN, pipenetwork::OPEN,
      pipenetwork::OPEN};

  std::vector<pipenetwork::Pipe_prop> pipe_props;
  for (int i = 0; i < pipe_ids.size(); ++i) {
    pipenetwork::Pipe_prop pipe_prop;
    pipe_prop.id = pipe_ids[i];
    pipe_prop.length = length[i];
    pipe_prop.diameter = diameter[i];
    pipe_prop.roughness = roughness[i];
    pipe_prop.node1_id = nodeids[i].first;
    pipe_prop.node2_id = nodeids[i].second;
    pipe_prop.status = status[i];

    pipe_props.emplace_back(pipe_prop);
  }

  mesh->create_pipes(pipe_props);

  //  mesh->print_summary ();

  double init_discharge = 1e-3;

  mesh->iterate_over_links(std::bind(&pipenetwork::Link::update_sim_discharge,
                                     std::placeholders::_1,
                                     init_discharge));  // initialze discharge
  bool pdd_mode = false;
  auto assembler = std::make_shared<pipenetwork::MatrixAssembler>(
      mesh, curves_info, pdd_mode);

  assembler->assemble_residual();
  assembler->update_jacobian();

  // A, x, b
  auto jac = assembler->jac_matrix();
  auto variables = assembler->variable_vector();
  auto residuals = assembler->residual_vector();

  SECTION("Pardiso 6 solver") {
    std::string solver_name = "pardiso";
    std::shared_ptr<pipenetwork::Solver> solve_ptr(
        Factory<pipenetwork::Solver>::instance()->create(solver_name));
    solve_ptr->assembled_matrices(jac, variables, residuals);
    auto x_diff = solve_ptr->solve();

    auto linear_sys_res = ((*jac) * x_diff - (*residuals)).norm();

    REQUIRE(linear_sys_res < tolerance);
  }

  SECTION("cuda solver") {
    std::string solver_name = "cuda";
    std::shared_ptr<pipenetwork::Solver> solve_ptr(
        Factory<pipenetwork::Solver>::instance()->create(solver_name));
    solve_ptr->assembled_matrices(jac, variables, residuals);
    auto x_diff = solve_ptr->solve();

    auto linear_sys_res = ((*jac) * x_diff - (*residuals)).norm();
    REQUIRE(linear_sys_res < tolerance);
  }
}
