#ifndef PIPE_NETWORK_PARDISO_UNSYM_H
#define PIPE_NETWORK_PARDISO_UNSYM_H

#include <exception>
#include <iostream>

#include <cmath>
#include <cstddef>

#include "solver.h"
#include <Eigen/Sparse>

namespace pipenetwork {

//! Pipe network Eigen GMRES class
//! \brief Eigen GMRES solver class using Eigen
    class Pardiso_unsym : public Solver {
    public:
        //! Constructor with tolerance, precision and iterations
        //! \param[in] max_iter Maximum number of iterations
        //! \param[in] tolerance Tolerance for solver to achieve convergence
        Pardiso_unsym(unsigned max_iter, double tolerance) : Solver(max_iter, tolerance) {

        }
        //! Call Pardiso Unsymmetric solver
        //! \retval status Return status of the solver
        bool solve() override;

    protected:
        //! Maximum number of iterations
        using Solver::max_iter_;
        //! Tolerance
        using Solver::tolerance_;
        //! Delta
        using Solver::delta_;
        //! Number of iterations
        using Solver::niterations_;
        //! Vector x
        using Solver::vec_x_;
        //! Vector b
        using Solver::vec_b_;
        //! Matrix A
        using Solver::mat_a_;


        //! Return the type of assembler for Factory
        //! \retval assembler_type Eigen assembler
        std::string assembler_type() const override { return "Pardiso_Unsym"; }
    };
}  // namespace pipenetwork

#endif //PIPE_NETWORK_PARDISO_UNSYM_H
