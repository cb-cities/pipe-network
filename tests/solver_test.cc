#include "catch.hpp"

#include "curves.h"
#include "factory.h"
//#include "input.h"
#include "matrix_assembler.h"
#include "solver.h"

// Check solver class
TEST_CASE("Solver is checked", "[Solver]") {

  // Tolerance
  const double tolerance = 1.e-8;

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

  std::vector<pipenetwork::JunctionProp> junc_props;
  for (int i = 0; i < elevations.size(); ++i) {
    pipenetwork::JunctionProp junc_prop;
    junc_prop.name = junction_ids[i];
    junc_prop.elevation = elevations[i];
    junc_prop.demand = demands[i];
    junc_prop.leak_diameter = leak_diameters[i];

    junc_props.emplace_back(junc_prop);
  }

  std::vector<std::string> res_ids{"13"};
  std::vector<double> heads{3.048e+02};

  std::vector<pipenetwork::ReservoirProp> res_props;
  for (int i = 0; i < res_ids.size(); ++i) {
    pipenetwork::ReservoirProp res_prop;
    res_prop.name = res_ids[i];
    res_prop.head = heads[i];

    res_props.emplace_back(res_prop);
  }

  std::vector<std::string> pipe_ids{"10", "11"};
  std::vector<std::pair<std::string, std::string>> nodeids{
      std::make_pair("13", "10"), std::make_pair("10", "11")};
  const std::vector<double> length{3209.5440000000003, 3209.5440000000003,
                                   1609.344, 1609.344};
  const std::vector<double> diameter{0.5588, 0.4572, 0.35559999999999997,
                                     0.254};
  const std::vector<double> roughness{100, 100, 100, 100};

  std::vector<pipenetwork::PipeProp> pipe_props;
  for (int i = 0; i < pipe_ids.size(); ++i) {
    pipenetwork::PipeProp pipe_prop;
    pipe_prop.name = pipe_ids[i];
    pipe_prop.length = length[i];
    pipe_prop.diameter = diameter[i];
    pipe_prop.roughness = roughness[i];
    pipe_prop.node1_name = nodeids[i].first;
    pipe_prop.node2_name = nodeids[i].second;

    pipe_props.emplace_back(pipe_prop);
  }

  std::vector<std::string> pump_ids{"12"};
  std::vector<std::pair<std::string, std::string>> pump_nodeids{
      std::make_pair("10", "12")};

  std::vector<pipenetwork::PumpProp> pump_props;

  for (int i = 0; i < pump_ids.size(); ++i) {
    pipenetwork::PumpProp pump_prop;
    pump_prop.name = pump_ids[i];
    pump_prop.type = pipenetwork::PumpType::POWERPUMP;
    pump_prop.node1_name = pump_nodeids[i].first;
    pump_prop.node2_name = pump_nodeids[i].second;
    pump_props.emplace_back(pump_prop);
  }

  std::vector<std::string> valve_ids{"13"};
  std::vector<std::pair<std::string, std::string>> valve_nodeids{
      std::make_pair("11", "12")};

  std::vector<pipenetwork::ValveProp> valve_props;

  for (int i = 0; i < valve_ids.size(); ++i) {
    pipenetwork::ValveProp valve_prop;
    valve_prop.name = valve_ids[i];
    valve_prop.type = pipenetwork::ValveType::PRVALVE;
    valve_prop.status = pipenetwork::LinkStatus ::ACTIVE;
    valve_prop.setting = 10;
    valve_prop.node1_name = valve_nodeids[i].first;
    valve_prop.node2_name = valve_nodeids[i].second;

    valve_props.emplace_back(valve_prop);
  }

  mesh->create_nodes(junc_props, res_props);
  mesh->create_links(pipe_props, pump_props, valve_props);
  mesh->create_mesh_graph();

  bool pdd_mode = false;
  auto assembler =
      std::make_shared<pipenetwork::linear_system::MatrixAssembler>(
          mesh, curves_info, pdd_mode);
  assembler->assemble_residual();
  assembler->update_jacobian();
  Eigen::VectorXd& var_vec = assembler->variable_vector();
  const Eigen::VectorXd& res_vec = assembler->residual_vector();
  const Eigen::SparseMatrix<double, Eigen::RowMajor>& jac_matrix =
      assembler->jac_matrix();

  SECTION("mkl pardiso solver") {
    std::string solver_name = "mkl_pardiso";
    std::shared_ptr<pipenetwork::linear_system::Solver> solve_ptr(
        Factory<pipenetwork::linear_system::Solver>::instance()->create(
            solver_name));

    solve_ptr->assembled_matrices(assembler);
    auto x_diff = solve_ptr->solve();
    auto linear_sys_res = (jac_matrix * x_diff - res_vec).norm();

    REQUIRE(linear_sys_res < tolerance);
  }

  SECTION("Eigen LU solver") {
    std::string solver_name = "LU";
    std::shared_ptr<pipenetwork::linear_system::Solver> solve_ptr(
        Factory<pipenetwork::linear_system::Solver>::instance()->create(
            solver_name));

    solve_ptr->assembled_matrices(assembler);
    auto x_diff = solve_ptr->solve();
    auto linear_sys_res = (jac_matrix * x_diff - res_vec).norm();

    REQUIRE(linear_sys_res < tolerance);
  }
}
