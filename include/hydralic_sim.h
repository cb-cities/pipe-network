#ifndef PIPE_NETWORK_HYDRALIC_SIM_H
#define PIPE_NETWORK_HYDRALIC_SIM_H
//#include <boost/filesystem.hpp>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#include "factory.h"
#include "io.h"
#include "matrix_assembler.h"
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
  //! \param[in] debug debug mode. if the mode is on, the initial variables,
  //! residual, and jacobian matrix will be recorded, iteration process will be
  //! printed
  Hydralic_sim(const std::shared_ptr<Mesh>& mesh,
               std::shared_ptr<Curves>& curves_info, bool pdd_mode = false,
               bool debug = false);

  //! run simulation
  //! \param[in] NR_tolerance residual tolerance for newton ralphson process
  //! \param[in] max_nr_steps max number of steps for newton ralphson process
  bool run_simulation(double NR_tolerance = 1.e-8, int max_nr_steps = 20);

  //! get the norm of simulation residual
  double sim_residual_norm() const {
    return assembler_->residual_vector().norm();
  }

  //! update mesh using simulation results
  void update_mesh();

 private:
  //! initial variable vector before linear system solving (for line search
  //! comparison)
  Eigen::VectorXd init_variable_;
  //! the mesh ptr
  std::shared_ptr<Mesh> mesh_;
  //! the assember ptr
  std::shared_ptr<linear_system::MatrixAssembler> assembler_;
  //! the solver ptr
  std::shared_ptr<linear_system::Solver> solver_;
  //! Network size Threshold for parallel solver
  Index nthre_{10000};
  //! debug flag
  bool debug_{false};
  //! backtracking (in line search) settings
  int bt_max_iter_{20};
  double bt_roh_{0.5};

  //! update system variable with line search
  void line_search_update(const Eigen::VectorXd& x_diff);
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_HYDRALIC_SIM_H
