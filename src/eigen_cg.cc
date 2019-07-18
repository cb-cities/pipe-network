#include "eigen_cg.h"

// Eigen CG Solver
bool pipenetwork::EigenCG::solve() {
    bool convergence = false;
    const size_t n = vec_b_->size();
    Eigen::VectorXd x_diff(n);

    // Use GMRES solver in Eigen to solve
    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double,Eigen::RowMajor> > lscg;
    lscg.setMaxIterations(max_iter_);
    lscg.setTolerance(tolerance_);

    lscg.compute(*mat_a_);
    x_diff = lscg.solve(*vec_b_);
    x_diff = x_diff.array();

    if (lscg.info() == 0) {
        convergence = true;
        //    std::cout << "difference :" <<x_diff<<std::endl;

        // Update x (variables) in the system
        (*vec_x_) -= x_diff;

        // Update delta and iterations
        this->delta_ = lscg.error();
        this->niterations_ = lscg.iterations();

        //    std::cout << "Iteration: " << niterations_ << " of " << max_iter_
        //              << "; delta: " << delta_ << std::endl;
    } else if (lscg.info() == 1) {
        throw std::runtime_error(
                "The provided data did not satisfy the prerequisites.");
    } else if (lscg.info() == 2) {
        std::cout << "Iteration: " << lscg.iterations() << " of " << max_iter_
                  << "; delta: " << lscg.error() << std::endl;
        throw std::runtime_error("Iterative procedure did not converge.");
    } else {
        throw std::runtime_error(
                "The inputs are invalid, or the algorithm has been improperly called.");
    }

    return convergence;
}
