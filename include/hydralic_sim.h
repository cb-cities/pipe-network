#ifndef PIPE_NETWORK_HYDRALIC_SIM_H
#define PIPE_NETWORK_HYDRALIC_SIM_H
#include <fstream>
#include <string>
#include <vector>

#include "factory.h"
#include "input.h"
#include "matrix_assembler.h"
#include "output.h"
#include "settings.h"
#include "solver.h"

namespace pipenetwork {
//! Hydraulic Simulation class
//! \brief Class for pipenetwork hydraulic simulation
class Hydralic_sim {
 public:
  //! Constructor with mesh input
  //! \param[in] mesh the mesh pointer
  //! \param[in] curves_info a pointer for curve information of the mesh
  //! \param[in] pdd_mode if simulation type is pressure demand driven or demand
  //! driven
  //! \param[in] solver_name name of the solver to use
  //! \param[in] debug debug mode. if the mode is on, the initial variables,
  //! residual, and jacobian matrix will be recorded, iteration process will be
  //! printed
  Hydralic_sim(const std::shared_ptr<Mesh>& mesh,
               std::shared_ptr<Curves>& curves_info, bool pdd_mode = false,
               const std::string& solver_name = "mkl_pardiso",
               bool debug = false) {
    mesh_ = mesh;
    assembler_ = std::make_shared<MatrixAssembler>(mesh, curves_info, pdd_mode);
    std::shared_ptr<Solver> solve_ptr(
        Factory<Solver>::instance()->create(solver_name));
    solver_ = solve_ptr;
    debug_ = debug;
  };

  //! Constructor with .inp file path
  //! \param[in] filepath .inp file path
  //! \param[in] mesh_name name of the mesh, used for creating result file name
  //! \param[in] pdd_mode if simulation type is pressure demand driven or demand
  //! driven
  //! \param[in] solver_name name of the solver to use
  //! \param[in] debug debug mode. if the mode is on, the initial variables,
  //! residual, and jacobian matrix will be recorded, iteration process will be
  //! printed
  Hydralic_sim(const std::string& filepath, const std::string& mesh_name,
               bool pdd_mode = false,
               const std::string& solver_name = "mkl_pardiso",
               bool debug = false);

  //! Constructor for synthetic net as input
  //! \param[in] syn_size the size of the synthesized network
  //! (syn_size,syn_size). \param[in] pdd_mode if simulation type is pressure
  //! demand driven or demand driven \param[in] solver_name name of the solver
  //! to use \param[in] debug debug mode. if the mode is on, the initial
  //! variables, residual, and jacobian matrix will be recorded, iteration
  //! process will be printed
  Hydralic_sim(int syn_size, bool pdd_mode,
               const std::string& solver_name = "mkl_pardiso",
               bool debug = false);

  //! run simulation
  //! \param[in] NR_tolerance residual tolerance for newton ralphson process
  //! \param[in] max_nr_steps max number of steps for newton ralphson process
  //! \param[in] line_search whether to use line search
  //! \param[in] output_path the path for output
  bool run_simulation(double NR_tolerance = 1.e-8, int max_nr_steps = 1000,
                      std::string output_path = "../results/res_",
                      bool line_search = true);
  //! get the norm of simulation residual
  double sim_residual_norm() const { return residual_norm_; }

  //! line search for newton ralphson
  //! \param[in] x_diff the variable difference vector from linear system solver
  void line_search_func(const Eigen::VectorXd& x_diff);

 private:
  //! original variable vector (for line search
  Eigen::VectorXd original_variable_;
  //! the mesh ptr
  std::shared_ptr<Mesh> mesh_;
  //! the assember ptr
  std::shared_ptr<MatrixAssembler> assembler_;
  //! the solver ptr
  std::shared_ptr<Solver> solver_;
  //! variable vector
  std::shared_ptr<Eigen::VectorXd> variables_;
  //! residual vector
  std::shared_ptr<Eigen::VectorXd> residuals_;
  //! solver tolerance
  double inner_solver_tolerance_{1e-8};
  //! iteration steps
  int max_solver_steps_{10000};
  //! residual norm
  double residual_norm_{std::numeric_limits<float>::min()};
  //! initial discharge
  double init_discharge_{1e-3};
  //! debug flag
  bool debug_{false};
  //! alpha for backtracking
  int bt_max_iter_{20};
  double bt_roh_{0.5};

  //! Function to write the final result
  //! \param[in] output_path the path for output
  //! \param[in] var variables need to be written
  void write_final_result(const std::string& output_path,
                          const Eigen::VectorXd& var);
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_HYDRALIC_SIM_H
