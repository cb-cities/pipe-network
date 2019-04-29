#include "catch.hpp"

#include "eigen_gmres.h"
#include "io.h"
#include "matrix_assembler.h"
#include "mesh.h"
#include "node.h"
#include "pipe.h"
#include "settings.h"

// Test to solve example network in Todini(2013)
TEST_CASE("Example network is tested", "[NR method]") {
    // Tolerance
    const double tolerance = 1.e-6;

    // Mesh index
    const unsigned meshid = 101;

    // Creat a mesh
    auto mesh = std::make_shared<pipenetwork::Mesh>(meshid);

    // Create a IO class object
    auto IO = std::make_unique<pipenetwork::IO>();

    // Read node and pipe information from CSV files
    bool read_network = IO->read_network("../benchmarks/todini_network_node.csv",
                                         "../benchmarks/todini_network_pipe.csv");
    std::vector<Eigen::Vector3d> coords = IO->nodal_coordinates();
    std::vector<std::pair<Index, double>> init_nodal_head =
            IO->initial_nodal_head();
    std::vector<std::pair<Index, double>> init_nodal_discharge =
            IO->initial_nodal_discharge();
    std::vector<std::pair<Index, Index>> node_pairs = IO->node_pairs();
    std::vector<double> diameter = IO->diameters();
    std::vector<double> roughness = IO->roughness();
    std::vector<bool> status = IO->pipe_status();
    std::vector<std::pair<Index, double>> init_pipe_discharge =
            IO->initial_pipe_discharge();

    double pipe_discharge_unknown = 30;
    // Create nodal pointers based on nodal coordinates in the mesh
    mesh->create_nodes(coords);

    // Create pipes based on pipe indices and previous created node pointers in
    // the mesh
    bool all_pipe_created =
            mesh->create_pipes(node_pairs, diameter, roughness, status);

    // Initialize pipe discharge
    mesh->initialize_pipe_discharge(pipe_discharge_unknown);

    // Assign initial nodal head and discharge
    mesh->assign_node_head(init_nodal_head);
    mesh->assign_node_discharge(init_nodal_discharge);

    // Initialize matrix assembler and obtain global index to nodes and pipes
    auto assembler = std::make_shared<pipenetwork::MatrixAssembler>();
    assembler->global_nodal_pipe_indices(mesh, 1);

    // Assemble variable vector, residual vector and Jacobian
    assembler->sim_assemble_variable_vector_v2();
    std::shared_ptr<Eigen::VectorXd> variable_vec = assembler->variable_vec();
    assembler->assemble_residual_vector_v2();
    std::shared_ptr<Eigen::VectorXd> residual_vec = assembler->residual_vec();
    assembler->sim_assemble_jacobian_v2();
    std::shared_ptr<Eigen::SparseMatrix<double>> jac = assembler->jac();

    // Creat a eigen gmres solver and solve
    double gmres_tolerance = 1.e-12;
    const unsigned max_iter = 3000;
    auto solver =
            std::make_shared<pipenetwork::EigenGMRES>(max_iter, gmres_tolerance);
    solver->assembled_matrices(jac, variable_vec, residual_vec);

    // Apply restraints
    Eigen::VectorXd restraints(11);
    restraints << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1;
    solver->restrains(restraints);

    std::cout << "Jac = " << std::endl
              << (*jac) << std::endl
              << std::endl
              << "residual = " << std::endl
              << std::endl
              << (*residual_vec) << std::endl
              << "variable = " << std::endl
              << (*variable_vec) << std::endl
              << std::endl;
    //
    bool issolved = solver->solve();
    assembler->sim_apply_variables_v2();
    assembler->assemble_residual_vector_v2();
    std::cout << "Jac = " << std::endl
              << (*jac) << std::endl
              << std::endl
              << "residual = " << std::endl
              << std::endl
              << (*residual_vec) << std::endl
              << "variable = " << std::endl
              << (*variable_vec) << std::endl
              << std::endl;

    unsigned nr_iter = 1;
    while (nr_iter < 20) {
        std::cout << "niter = " << nr_iter << std::endl;
        bool issolved = solver->solve();
        assembler->sim_apply_variables_v2();
        assembler->assemble_residual_vector_v2();
        //    assembler->sim_assemble_variable_vector_v2();
        assembler->sim_assemble_jacobian_v2();
        if (residual_vec->norm() < tolerance) {
            std::cout << "Jac = "
                      << std::endl
                      //                << (*jac) << std::endl
                      //                << std::endl
                      << "residual = " << std::endl
                      << (*residual_vec) << std::endl
                      << std::endl
                      << "variable = " << std::endl
                      << (*variable_vec) << std::endl
                      << std::endl;
            break;
        }
        std::cout << "Jac = "
                  << std::endl
                  //              << (*jac) << std::endl
                  //              << std::endl
                  << "residual = " << std::endl
                  << (*residual_vec) << std::endl
                  << std::endl
                  << "variable = " << std::endl
                  << (*variable_vec) << std::endl
                  << std::endl;
        //
        nr_iter++;
    }
    REQUIRE(residual_vec->norm() == Approx(0.0).epsilon(tolerance));
}
