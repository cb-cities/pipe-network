#ifndef PIPE_NETWORK_HYDRALIC_SIM_H
#define PIPE_NETWORK_HYDRALIC_SIM_H
#include <fstream>
#include <string>
#include <vector>

#include "eigen_cg.h"
#include "eigen_gmres.h"
#include "input.h"
#include "matrix_assembler.h"
#include "settings.h"

namespace pipenetwork {
//! Hydraulic Simulation class
//! \brief Class for pipenetwork hydraulic simulation
class Hydralic_sim {
 public:
  explicit Hydralic_sim(const std::shared_ptr<Mesh>& mesh,
                        bool pdd_mode = false, bool debug = false) {
    assembler_ = std::make_shared<MatrixAssembler>(mesh, pdd_mode);
    solver_ =
        std::make_shared<EigenCG>(max_solver_steps_, inner_solver_tolerance_);
    debug_ = debug;
  };
  Hydralic_sim(const std::string& filepath,
               const std::vector<double>& leak_diameters, bool pdd_mode = false,
               bool debug = false);

  //! run simulation
  bool run_simulation(double NR_tolerance = 1.e-8, int max_nr_steps = 100);
  //! get the norm of simulation residual
  double sim_residual_norm() const { return residual_norm_; }

 private:
  //! the assember ptr
  std::shared_ptr<MatrixAssembler> assembler_;
  //! the solver ptr
  std::shared_ptr<EigenCG> solver_;
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
};
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_HYDRALIC_SIM_H
