#ifndef PIPE_NETWORK_HYDRALIC_SIM_H
#define PIPE_NETWORK_HYDRALIC_SIM_H

#include "eigen_gmres.h"
#include "matrix_assembler.h"
#include "input.h"
#include "settings.h"

namespace pipenetwork {
//! Hydraulic Simulation class
//! \brief Class for pipenetwork hydraulic simulation
class Hydralic_sim {
 public:
  explicit Hydralic_sim(const std::shared_ptr<Mesh>& mesh, bool pdd_mode = false) {
    assembler_ = std::make_shared<MatrixAssembler>(mesh, pdd_mode);
    solver_ = std::make_shared<EigenGMRES>(max_solver_steps_,
                                           inner_solver_tolerance_);
  };
  Hydralic_sim(const std::string & filepath, const std::vector<double> & leak_diameters, bool pdd_mode = false);

  //! run simulation
  bool run_simulation(double NR_tolerance = 1.e-8, int max_nr_steps = 100);
  //! get the norm of simulation residual
  double sim_residual_norm() const { return residual_norm_; }

 private:
  //! the assember ptr
  std::shared_ptr<MatrixAssembler> assembler_;
  //! the solver ptr
  std::shared_ptr<EigenGMRES> solver_;
  //! solver tolerance
  double inner_solver_tolerance_{1e-12};
  //! iteration steps
  int max_solver_steps_{5000};
  //! residual norm
  double residual_norm_{-9999};
  //! initial discharge
  double init_discharge_{1e-3};
};
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_HYDRALIC_SIM_H
