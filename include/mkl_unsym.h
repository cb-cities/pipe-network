#ifndef PIPE_NETWORK_MKL_UNSYM_H
#define PIPE_NETWORK_MKL_UNSYM_H

#include <exception>
#include <iostream>
#include <cmath>
#include <cstddef>
#include <Eigen/Sparse>

#include "solver.h"
#include "mkl_pardiso.h"


namespace pipenetwork {

//! Pipe network Eigen GMRES class
//! \brief Eigen GMRES solver class using Eigen
    class Mkl_unsym : public Solver {
    public:
        //! Constructor with tolerance, precision and iterations
        //! \param[in] max_iter Maximum number of iterations
        //! \param[in] tolerance Tolerance for solver to achieve convergence
        Mkl_unsym();
        //! Call Pardiso Unsymmetric solver
        //! \retval status Return status of the solver
        Eigen::VectorXd solve() override;

    protected:
        //! Vector x
        using Solver::vec_x_;
        //! Vector b
        using Solver::vec_b_;
        //! Matrix A
        using Solver::mat_a_;
        //! Row index for csr format matrix
        using Solver::ia_;
        //! Column index for csr format matrix
        using Solver::ja_;
        //! Matrix values for csr format matrix
        using Solver::a_;
        //! Return the type of assembler for Factory
        //! \retval assembler_type Eigen assembler
        std::string assembler_type() const override { return "Pardiso_Unsym"; }

    private:
        int mtype_ = 11; /* Real unsymmetric matrix */
        /* RHS and solution vectors. */
        int nrhs_ = 1; /* Number of right hand sides. */

        /* Internal solver memory pointer pt,                  */
        /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
        /* or void *pt[64] should be OK on both architectures  */
        void* pt_[64];
        /* Pardiso control parameters. */
        int iparm_[64];
        int maxfct_, mnum_, phase_, error_, msglvl_;
        /* Number of processors. */
        int num_procs_;
        double ddum_; /* Double dummy */
        int idum_;    /* Integer dummy. */
    };
}  // namespace pipenetwork

#endif //PIPE_NETWORK_MKL_UNSYM_H
