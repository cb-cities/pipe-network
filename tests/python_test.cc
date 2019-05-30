#include "catch.hpp"

#include "eigen_gmres.h"
#include "input.h"
#include "matrix_assembler.h"
#include "mesh.h"
#include "node.h"
#include "pipe.h"
#include "settings.h"


// Test to solve example network in Todini(2013)
TEST_CASE("Python example network is tested", "[NR method]") {
    // Tolerance
    const double tolerance = 1.e-6;

    // Mesh index
    const unsigned meshid = 101;

    // Creat a mesh
    auto mesh = std::make_shared<pipenetwork::Mesh>(meshid);

    // Read node and pipe information from INP files
    auto IO = std::make_unique<pipenetwork::Input>("../benchmarks/Net1c.inp");

    auto coords = IO->get_node_coord ();
    auto node_elevation = IO->get_node_elevation ();
    auto node_demand = IO->get_node_demand ();

    auto pipe_end_nodes = IO->get_pipe_end_ids ();
    auto rough = IO->get_pipe_roughness ();
    auto dia = IO->get_pipe_diameters ();
    auto status = IO->get_pipe_status ();
    auto length = IO->get_pipe_length ();



    double pipe_discharge_unknown = 1e-3;
    // Create nodal pointers based on nodal coordinates in the mesh
    mesh->create_nodes(coords);

    // Create pipes based on pipe indices and previous created node pointers in
    // the mesh
    bool all_pipe_created =
            mesh->create_pipes(pipe_end_nodes, dia,length, rough, status);

    // Initialize pipe discharge
    mesh->initialize_pipe_discharge(pipe_discharge_unknown);

    // Assign initial nodal head and discharge
    mesh->assign_node_elevation(node_elevation);
    mesh->assign_node_demand(node_demand);


    std::cout<<junction<<pipe<<std::endl;




//    double leak_dia1 = 0.5;
//    double leak_dia2 = 0.2;
//    std::vector<std::pair<Index, double>> node_leak_dia;
//    node_leak_dia.emplace_back(std::make_pair(1, leak_dia1));
//    node_leak_dia.emplace_back(std::make_pair(2, leak_dia2));
//
//    mesh->assign_node_leak (node_leak_dia);
//
//    SECTION("DD mode test") {
//        // Initialize matrix assembler and obtain global index to nodes and pipes
//        bool pdd_mode = false;
//        auto assembler = std::make_shared<pipenetwork::MatrixAssembler>(pdd_mode);
//        assembler->global_nodal_pipe_indices(mesh);
//
//        // Assemble variable vector, residual vector and Jacobian
//        assembler->assemble_variable_vector();
//        std::shared_ptr<Eigen::VectorXd> variable_vec = assembler->variable_vec();
//        assembler->assemble_residual_vector();
//        std::shared_ptr<Eigen::VectorXd> residual_vec = assembler->residual_vec();
//        assembler->assemble_jacobian();
//        std::shared_ptr<Eigen::SparseMatrix<double>> jac = assembler->jac();
//
//        // Creat a eigen gmres solver and solve
//        double gmres_tolerance = 1.e-12;
//        const unsigned max_iter = 10000;
//        auto solver =
//                std::make_shared<pipenetwork::EigenGMRES>(max_iter, gmres_tolerance);
//        solver->assembled_matrices(jac, variable_vec, residual_vec);
//
//
//        unsigned niter = 10;
//        for (unsigned nr_iter = 0; nr_iter < niter; ++nr_iter) {
////            std::cout << "niter = " << nr_iter << std::endl;
////            std::cout << "Jac = " << std::endl
////                      << (*jac) << std::endl
////                      << std::endl
////                      << "residual = " << std::endl
////                      << (*residual_vec) << std::endl
////                      << std::endl
////                      << "variable = " << std::endl
////                      << (*variable_vec) << std::endl
////                      << std::endl;
//
//            bool issolved = solver->solve();
//            assembler->apply_variables();
//            assembler->assemble_residual_vector();
//            assembler->assemble_jacobian();
//            std::cout << "niter = " << nr_iter << std::endl;
//            std::cout<< "residual = " << std::endl
//                      << residual_vec->norm() << std::endl;
//
//
//            if (residual_vec->norm() < tolerance) {
////                std::cout << "niter = " << nr_iter << std::endl;
////                std::cout << "Jac = " << std::endl
////                          << (*jac) << std::endl
////                          << std::endl
////                          << "residual = " << std::endl
////                          << (*residual_vec) << std::endl
////                          << std::endl
////                          << "variable = " << std::endl
////                          << (*variable_vec) << std::endl
////                          << std::endl;
////                break;
//            }
//        }
////        std::cout<<(*variable_vec) << std::endl;
////        REQUIRE(residual_vec->norm() == Approx(0.0).epsilon(tolerance));
//    }

//    SECTION("PD mode test") {
//        // Initialize matrix assembler and obtain global index to nodes and pipes
//        bool pdd_mode = true;
//        auto assembler = std::make_shared<pipenetwork::MatrixAssembler>(pdd_mode);
//        assembler->global_nodal_pipe_indices(mesh);
//
//        // Assemble variable vector, residual vector and Jacobian
//        assembler->assemble_variable_vector();
//        std::shared_ptr<Eigen::VectorXd> variable_vec = assembler->variable_vec();
//        assembler->assemble_residual_vector();
//        std::shared_ptr<Eigen::VectorXd> residual_vec = assembler->residual_vec();
//        assembler->assemble_jacobian();
//        std::shared_ptr<Eigen::SparseMatrix<double>> jac = assembler->jac();
//
//        // Creat a eigen gmres solver and solve
//        double gmres_tolerance = 1.e-12;
//        const unsigned max_iter = 5000;
//        auto solver =
//                std::make_shared<pipenetwork::EigenGMRES>(max_iter, gmres_tolerance);
//        solver->assembled_matrices(jac, variable_vec, residual_vec);
//
//        bool issolved = solver->solve();
//        assembler->apply_variables();
//        assembler->assemble_residual_vector();
//
//        unsigned niter = 30;
//        for (unsigned nr_iter = 0; nr_iter < niter; ++nr_iter) {
//
//            bool issolved = solver->solve();
//            assembler->apply_variables();
//            assembler->assemble_residual_vector();
//            assembler->assemble_jacobian();
//            if (residual_vec->norm() < tolerance) {
//                std::cout << "niter = " << nr_iter << std::endl;
//                std::cout << "Jac = " << std::endl
//                          << (*jac) << std::endl
//                          << std::endl
//                          << "residual = " << std::endl
//                          << (*residual_vec) << std::endl
//                          << std::endl
//                          << "variable = " << std::endl
//                          << (*variable_vec) << std::endl
//                          << std::endl;
//                break;
//            }
//        }
//        REQUIRE(residual_vec->norm() == Approx(0.0).epsilon(tolerance));
//    }
}
