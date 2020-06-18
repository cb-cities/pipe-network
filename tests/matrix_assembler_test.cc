#include "catch.hpp"
#include <cmath>
#include <iomanip>

#include "curves.h"
//#include "input.h"
#include "matrix_assembler.h"

// Check matrix_assembler class
TEST_CASE("MatrixAssembler is checked", "[MatrixAssembler]") {

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

  //  mesh->print_summary();

  SECTION("TEST VARIABLES") {
    auto variables = pipenetwork::linear_system::Variables(mesh);
    SECTION("variable vector") {
      auto var_vec = variables.variables_vec();
      REQUIRE(var_vec.size() == 13);
      REQUIRE(var_vec[0] == Approx(216.408).epsilon(tolerance));
      REQUIRE(var_vec[5] == Approx(9.464e-03).epsilon(tolerance));
      REQUIRE(var_vec[10] ==
              Approx(pipenetwork::INIT_FLOWRATE).epsilon(tolerance));
      REQUIRE(var_vec[12] == Approx(0).epsilon(tolerance));
    }
    SECTION("variable vector") {
      auto demands_heads_vec = variables.demands_heads_vec();
      REQUIRE(demands_heads_vec.size() == 4);
      REQUIRE(demands_heads_vec[0] == Approx(0).epsilon(tolerance));
      REQUIRE(demands_heads_vec[1] == Approx(9.464e-03).epsilon(tolerance));
      REQUIRE(demands_heads_vec[2] == Approx(9.464e-03).epsilon(tolerance));
      REQUIRE(demands_heads_vec[3] == Approx(3.048e+02).epsilon(tolerance));
    }
    SECTION("elevation vector") {
      auto elevations_vec = variables.elevations();
      REQUIRE(elevations_vec.size() == 4);
      REQUIRE(elevations_vec[0] == Approx(216.408).epsilon(tolerance));
      REQUIRE(elevations_vec[1] == Approx(216.408).epsilon(tolerance));
      REQUIRE(elevations_vec[2] == Approx(213.36).epsilon(tolerance));
      REQUIRE(elevations_vec[3] == Approx(3.048e+02).epsilon(tolerance));
    }
    SECTION("resistance vector") {
      auto link_resistance_coeff_vec = variables.link_resistance_coeff_vec();

      double coeff1 = pipenetwork::HW_COEFF * std::pow(roughness[1], -1.852) *
                      std::pow(diameter[1], -4.871) * length[1];

      REQUIRE(link_resistance_coeff_vec.size() == 4);
      REQUIRE(link_resistance_coeff_vec[1] ==
              Approx(coeff1).epsilon(tolerance));
    }

    SECTION("minor vector") {
      auto link_minor_loss_coeff_vec = variables.link_minor_loss_coeff_vec();
      REQUIRE(link_minor_loss_coeff_vec.size() == 4);
      REQUIRE(link_minor_loss_coeff_vec[0] == Approx(0).epsilon(tolerance));
      REQUIRE(link_minor_loss_coeff_vec[1] == Approx(0).epsilon(tolerance));
      REQUIRE(link_minor_loss_coeff_vec[2] == Approx(0).epsilon(tolerance));
      REQUIRE(link_minor_loss_coeff_vec[3] == Approx(0).epsilon(tolerance));
    }

    SECTION("leak area vector") {
      auto leak_areas = variables.leak_areas();
      double leak_area = pipenetwork::PI * std::pow((0.1 / 2), 2);
      REQUIRE(leak_areas.size() == 1);
      REQUIRE(leak_areas[0] == Approx(leak_area).epsilon(tolerance));
    }
  }

  SECTION("TEST RESIDUALS") {
    auto variables =
        std::make_shared<pipenetwork::linear_system::Variables>(mesh);
    auto residuals =
        pipenetwork::linear_system::Residuals(mesh, variables, curves_info);
    residuals.assemble_residual();
    auto res_vec = residuals.residual_vec();

    SECTION("node balance residual") {
      REQUIRE(res_vec.coeff(0) == Approx(-0.001).epsilon(tolerance));
      REQUIRE(res_vec.coeff(1) == Approx(-9.464e-03).epsilon(tolerance));
    }

    SECTION("node demand head residual") {
      REQUIRE(res_vec.coeff(5) == Approx(0).epsilon(tolerance));
      REQUIRE(res_vec.coeff(7) == Approx(0).epsilon(tolerance));
    }

    SECTION("pipe headloss residual") {
      // pipe headloss
      REQUIRE(res_vec.coeff(8) == Approx(-88.3916796736).epsilon(tolerance));
      // pump headloss
      double pump_res = 50 + (elevations[0] - elevations[2]) *
                                 pipenetwork::INIT_FLOWRATE * pipenetwork::G *
                                 1000.0;
      REQUIRE(res_vec.coeff(10) == Approx(pump_res).epsilon(tolerance));

      // valve headloss
      REQUIRE(res_vec.coeff(11) == Approx(-10).epsilon(tolerance));
    }

    SECTION("pipe leak residual") {
      REQUIRE(res_vec.coeff(12) == Approx(0).epsilon(tolerance));
    }
  }
  SECTION("TEST PDD RESIDUALS") {
    auto variables =
        std::make_shared<pipenetwork::linear_system::Variables>(mesh);
    auto residuals =
        pipenetwork::linear_system::Residuals(mesh, variables, curves_info);
    residuals.assemble_residual_pdd();
    auto res_vec = residuals.residual_vec();

    SECTION("node demand head residual") {
      REQUIRE(res_vec.coeff(5) == Approx(9.464e-03).epsilon(tolerance));
      REQUIRE(res_vec.coeff(7) == Approx(0).epsilon(tolerance));
    }
  }

  SECTION("TEST Jacobian") {
    auto variables =
        std::make_shared<pipenetwork::linear_system::Variables>(mesh);
    auto residuals = std::make_shared<pipenetwork::linear_system::Residuals>(
        mesh, variables, curves_info);
    residuals->assemble_residual();
    auto jacobian = pipenetwork::linear_system::Jacobian(
        mesh, variables, residuals, curves_info);

    SECTION("jacobian test: DD mode") {
      jacobian.update_jacobian();
      auto jac_matrix = jacobian.jac_matrix();
      // sub jacobian A: node_bal equation with respect to demand
      REQUIRE(jac_matrix.coeff(0, 4) == -1);
      REQUIRE(jac_matrix.coeff(1, 5) == -1);
      REQUIRE(jac_matrix.coeff(2, 6) == -1);
      REQUIRE(jac_matrix.coeff(3, 7) == -1);
      // sub jacobian B: node_bal equation with respect to flow
      REQUIRE(jac_matrix.coeff(3, 8) == -1);
      REQUIRE(jac_matrix.coeff(0, 8) == 1);
      REQUIRE(jac_matrix.coeff(0, 9) == -1);
      REQUIRE(jac_matrix.coeff(1, 9) == 1);
      REQUIRE(jac_matrix.coeff(0, 10) == -1);
      REQUIRE(jac_matrix.coeff(2, 10) == 1);
      REQUIRE(jac_matrix.coeff(1, 11) == -1);
      REQUIRE(jac_matrix.coeff(2, 11) == 1);
      // sub jacobian C: node_bal equation with respect to  leak flow
      REQUIRE(jac_matrix.coeff(0, 12) == 0);
      REQUIRE(jac_matrix.coeff(1, 12) == -1);
      REQUIRE(jac_matrix.coeff(2, 12) == 0);
      REQUIRE(jac_matrix.coeff(3, 12) == 0);
      // sub jacobian D: demand/head equation with respect to head
      REQUIRE(jac_matrix.coeff(4, 0) == 0);
      REQUIRE(jac_matrix.coeff(5, 1) == 0);
      REQUIRE(jac_matrix.coeff(6, 2) == 0);
      REQUIRE(jac_matrix.coeff(7, 3) == 1);
      // sub jacobian E: demand/head equation with respect to flow
      REQUIRE(jac_matrix.coeff(4, 4) == 1);
      REQUIRE(jac_matrix.coeff(5, 5) == 1);
      REQUIRE(jac_matrix.coeff(6, 6) == 1);
      REQUIRE(jac_matrix.coeff(7, 7) == 0);
      // sub jacobian F: headloss equation with respect to head
      REQUIRE(jac_matrix.coeff(8, 3) == -1);
      REQUIRE(jac_matrix.coeff(8, 0) == 1);
      REQUIRE(jac_matrix.coeff(9, 0) == -1);
      REQUIRE(jac_matrix.coeff(9, 1) == 1);
      REQUIRE(jac_matrix.coeff(10, 0) ==
              1000.0 * pipenetwork::G * pipenetwork::INIT_FLOWRATE);
      REQUIRE(jac_matrix.coeff(10, 2) ==
              -1000.0 * pipenetwork::G * pipenetwork::INIT_FLOWRATE);
      REQUIRE(jac_matrix.coeff(11, 1) == 0);
      REQUIRE(jac_matrix.coeff(11, 2) == 1);
      // sub jacobian G: headloss equation with respect to flow
      REQUIRE(jac_matrix.coeff(8, 8) ==
              Approx(0.5932445506024201).epsilon(tolerance));
      REQUIRE(jac_matrix.coeff(9, 9) ==
              Approx(1.576675297968889).epsilon(tolerance));
      auto pump_jac = (elevations[0] - elevations[2]) * 1000 * pipenetwork::G;
      REQUIRE(jac_matrix.coeff(10, 10) == Approx(pump_jac).epsilon(tolerance));
      REQUIRE(jac_matrix.coeff(11, 11) == 0);  // PRV
      // sub jacobian H: leak flow to flow
      REQUIRE(jac_matrix.coeff(12, 0) == 0);
      REQUIRE(jac_matrix.coeff(12, 1) == -1e-4);
      // sub jacobian I: leak flow to leak flow
      REQUIRE(jac_matrix.coeff(12, 12) == 1);
    }
    SECTION("jacobian test: PDD mode") {
      jacobian.update_jacobian_pdd();
      auto jac_matrix = jacobian.jac_matrix();

      // sub jacobian D: demand/head equation with respect to head
      REQUIRE(jac_matrix.coeff(4, 0) == 0);
      REQUIRE(jac_matrix.coeff(5, 1) == Approx(0).epsilon(tolerance));
      REQUIRE(jac_matrix.coeff(6, 2) == Approx(0).epsilon(tolerance));
      REQUIRE(jac_matrix.coeff(7, 3) == 1);
    }
  }

  SECTION("TEST Matrix Assembler") {
    bool pdd_mode = false;
    auto assembler =
        std::make_shared<pipenetwork::linear_system::MatrixAssembler>(
            mesh, curves_info, pdd_mode);
    assembler->system_update();
    Eigen::VectorXd& var_vec = assembler->variable_vector();
    auto res_vec = assembler->residual_vector();
    auto jac_matrix = assembler->jac_matrix();
    SECTION("variable vector") {
      REQUIRE(var_vec.size() == 13);
      REQUIRE(var_vec[0] == Approx(216.408).epsilon(tolerance));
      REQUIRE(var_vec[5] == Approx(9.464e-03).epsilon(tolerance));
      REQUIRE(var_vec[10] ==
              Approx(pipenetwork::INIT_FLOWRATE).epsilon(tolerance));
      REQUIRE(var_vec[12] == Approx(0).epsilon(tolerance));
    }
    SECTION("node balance residual") {
      REQUIRE(res_vec.coeff(0) == Approx(-0.001).epsilon(tolerance));
      REQUIRE(res_vec.coeff(1) == Approx(-9.464e-03).epsilon(tolerance));
      REQUIRE(res_vec.coeff(5) == Approx(0).epsilon(tolerance));
      REQUIRE(res_vec.coeff(7) == Approx(0).epsilon(tolerance));
      REQUIRE(res_vec.coeff(8) == Approx(-88.3916796736).epsilon(tolerance));
      // pump headloss
      double pump_res = 50 + (elevations[0] - elevations[2]) *
                                 pipenetwork::INIT_FLOWRATE * pipenetwork::G *
                                 1000.0;
      REQUIRE(res_vec.coeff(10) == Approx(pump_res).epsilon(tolerance));

      // valve headloss
      REQUIRE(res_vec.coeff(11) == Approx(-10).epsilon(tolerance));
    }
    SECTION("jacobian test: DD mode") {
      REQUIRE(jac_matrix.coeff(0, 4) == -1);
      REQUIRE(jac_matrix.coeff(1, 5) == -1);
      // sub jacobian A: node_bal equation with respect to demand
      REQUIRE(jac_matrix.coeff(0, 4) == -1);
      REQUIRE(jac_matrix.coeff(1, 5) == -1);
      REQUIRE(jac_matrix.coeff(2, 6) == -1);
      REQUIRE(jac_matrix.coeff(3, 7) == -1);
      // sub jacobian B: node_bal equation with respect to flow
      REQUIRE(jac_matrix.coeff(3, 8) == -1);
      REQUIRE(jac_matrix.coeff(0, 8) == 1);
      REQUIRE(jac_matrix.coeff(0, 9) == -1);
      REQUIRE(jac_matrix.coeff(1, 9) == 1);
      REQUIRE(jac_matrix.coeff(0, 10) == -1);
      REQUIRE(jac_matrix.coeff(2, 10) == 1);
      REQUIRE(jac_matrix.coeff(1, 11) == -1);
      REQUIRE(jac_matrix.coeff(2, 11) == 1);
      // sub jacobian C: node_bal equation with respect to  leak flow
      REQUIRE(jac_matrix.coeff(0, 12) == 0);
      REQUIRE(jac_matrix.coeff(1, 12) == -1);
      REQUIRE(jac_matrix.coeff(2, 12) == 0);
      REQUIRE(jac_matrix.coeff(3, 12) == 0);
      // sub jacobian D: demand/head equation with respect to head
      REQUIRE(jac_matrix.coeff(4, 0) == 0);
      REQUIRE(jac_matrix.coeff(5, 1) == 0);
      REQUIRE(jac_matrix.coeff(6, 2) == 0);
      REQUIRE(jac_matrix.coeff(7, 3) == 1);
      // sub jacobian E: demand/head equation with respect to flow
      REQUIRE(jac_matrix.coeff(4, 4) == 1);
      REQUIRE(jac_matrix.coeff(5, 5) == 1);
      REQUIRE(jac_matrix.coeff(6, 6) == 1);
      REQUIRE(jac_matrix.coeff(7, 7) == 0);
      // sub jacobian F: headloss equation with respect to head
      REQUIRE(jac_matrix.coeff(8, 3) == -1);
      REQUIRE(jac_matrix.coeff(8, 0) == 1);
      REQUIRE(jac_matrix.coeff(9, 0) == -1);
      REQUIRE(jac_matrix.coeff(9, 1) == 1);
      REQUIRE(jac_matrix.coeff(10, 0) ==
              1000.0 * pipenetwork::G * pipenetwork::INIT_FLOWRATE);
      REQUIRE(jac_matrix.coeff(10, 2) ==
              -1000.0 * pipenetwork::G * pipenetwork::INIT_FLOWRATE);
      REQUIRE(jac_matrix.coeff(11, 1) == 0);
      REQUIRE(jac_matrix.coeff(11, 2) == 1);
      // sub jacobian G: headloss equation with respect to flow
      REQUIRE(jac_matrix.coeff(8, 8) ==
              Approx(0.5932445506024201).epsilon(tolerance));
      REQUIRE(jac_matrix.coeff(9, 9) ==
              Approx(1.576675297968889).epsilon(tolerance));
      auto pump_jac = (elevations[0] - elevations[2]) * 1000 * pipenetwork::G;
      REQUIRE(jac_matrix.coeff(10, 10) == Approx(pump_jac).epsilon(tolerance));
      REQUIRE(jac_matrix.coeff(11, 11) == 0);  // PRV
      // sub jacobian H: leak flow to flow
      REQUIRE(jac_matrix.coeff(12, 0) == 0);
      REQUIRE(jac_matrix.coeff(12, 1) == -1e-4);
      // sub jacobian I: leak flow to leak flow
      REQUIRE(jac_matrix.coeff(12, 12) == 1);
    }
    SECTION("Test variable update") {
      Eigen::VectorXd& var_vec_test = assembler->variable_vector();
      Eigen::VectorXd x_diff(var_vec_test.size());

      x_diff.setOnes();
      var_vec_test -= x_diff;
      auto updated_var_vec = assembler->variable_vector();
      REQUIRE(updated_var_vec[0] == Approx(216.408 - 1).epsilon(tolerance));
      REQUIRE(updated_var_vec[5] == Approx(9.464e-03 - 1).epsilon(tolerance));
      REQUIRE(updated_var_vec[10] ==
              Approx(pipenetwork::INIT_FLOWRATE - 1).epsilon(tolerance));
      REQUIRE(updated_var_vec[12] == Approx(0 - 1).epsilon(tolerance));

      assembler->assemble_residual();
      auto updated_res_vec = assembler->residual_vector();
      REQUIRE(updated_res_vec.coeff(5) == Approx(-1).epsilon(tolerance));
      REQUIRE(updated_res_vec.coeff(7) == Approx(-1).epsilon(tolerance));

      assembler->update_jacobian();
      auto updated_jac_matrix = assembler->jac_matrix();
      // sub jacobian G: headloss equation with respect to flow
      REQUIRE(updated_jac_matrix.coeff(8, 8) ==
              Approx(213.2374859317).epsilon(tolerance));
      REQUIRE(updated_jac_matrix.coeff(11, 11) == 0);  // PRV
      // sub jacobian H: leak flow to flow
      REQUIRE(updated_jac_matrix.coeff(12, 0) == 0);
      REQUIRE(updated_jac_matrix.coeff(12, 1) == -1e-4);

      //          REQUIRE(updated_res_vec.coeff(8) ==
      //          Approx(-88.3916796736).epsilon(tolerance));
    }
  }
}