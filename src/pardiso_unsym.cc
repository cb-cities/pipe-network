#include "pardiso_unsym.h"
#include "factory.h"

using namespace std;
static Register<pipenetwork::Solver, pipenetwork::Pardiso_unsym> registry(
    "pardiso");

/* PARDISO prototype. */
extern "C" void pardisoinit(void*, int*, int*, int*, double*, int*);
extern "C" void pardiso(void*, int*, int*, int*, int*, int*, double*, int*,
                        int*, int*, int*, int*, int*, double*, double*, int*,
                        double*);
extern "C" void pardiso_chkmatrix(int*, int*, double*, int*, int*, int*);
extern "C" void pardiso_chkvec(int*, int*, double*, int*);
extern "C" void pardiso_printstats(int*, int*, double*, int*, int*, int*,
                                   double*, int*);

pipenetwork::Pardiso_unsym::Pardiso_unsym() : Solver() {
  // configure pardiso
  /* Auxiliary variables. */
  char* var;

  mtype_ = 11; /* Real unsymmetric matrix */
  /* RHS and solution vectors. */
  nrhs_ = 1; /* Number of right hand sides. */
  /* -------------------------------------------------------------------- */
  /* ..  Setup Pardiso control parameters and initialize the solvers      */
  /*     internal adress pointers. This is only necessary for the FIRST   */
  /*     call of the PARDISO solver.                                      */
  /* ---------------------------------------------------------------------*/
  error_ = 0;
  solver_ = 0; /* use sparse direct solver */
  pardisoinit(pt_, &mtype_, &solver_, iparm_, dparm_, &error_);

  if (error_ != 0) {
    if (error_ == -10) printf("No license file found \n");
    if (error_ == -11) printf("License is expired \n");
    if (error_ == -12) printf("Wrong username or hostname \n");
    exit(1);
  } else
    printf("[PARDISO]: License check was successful ... \n");

  /* Numbers of processors, value of OMP_NUM_THREADS */
  var = getenv("OMP_NUM_THREADS");
  if (var != NULL)
    sscanf(var, "%d", &num_procs_);
  else {
    num_procs_ = omp_get_max_threads();
  }

  /* -------------------------------------------------------------------- */
  /* .. Setup Pardiso control parameters. */
  /* -------------------------------------------------------------------- */
  for (int i = 0; i < 64; i++) {
    iparm_[i] = 0;
  }
  iparm_[0] = 1; /* No solver default */
  iparm_[1] = 1; /* Fill-in reordering from METIS */
  /* Numbers of processors, value of OMP_NUM_THREADS */
  iparm_[2] = num_procs_;
  iparm_[3] = 0;   /* No iterative-direct algorithm */
  iparm_[4] = 0;   /* No user fill-in reducing permutation */
  iparm_[5] = 0;   /* Write solution into x */
  iparm_[6] = 0;   /* Not in use */
  iparm_[7] = 2;   /* Max numbers of iterative refinement steps */
  iparm_[8] = 0;   /* Not in use */
  iparm_[9] = 13;   /* Perturb the pivot elements with 1E-13 */
  iparm_[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
  iparm_[11] = 0;  /* Not in use */
  iparm_[12] = 0;  /* Not in use */
  iparm_[13] = 0;  /* Output: Number of perturbed pivots */
  iparm_[14] = 0;  /* Not in use */
  iparm_[15] = 0;  /* Not in use */
  iparm_[16] = 0;  /* Not in use */
  iparm_[17] = -1; /* Output: Number of nonzeros in the factor LU */
  iparm_[18] = -1; /* Output: Mflops for LU factorization */
  iparm_[19] = 0;  /* Output: Numbers of CG Iterations */
  maxfct_ = 1;     /* Maximum number of numerical factorizations. */
  mnum_ = 1;       /* Which factorization to use. */
  msglvl_ = 0;     /* Print statistical information in file */
  error_ = 0;      /* Initialize error flag */
}

Eigen::VectorXd pipenetwork::Pardiso_unsym::solve() {
  // configure matrix
  int n = vec_b_->size();
  double* vec_b = vec_b_->data();
  double x_diff[n];
  int nnz = ia_[n];

  /* -------------------------------------------------------------------- */
  /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
  /*     notation.                                                        */
  /* -------------------------------------------------------------------- */
  for (int i = 0; i < n + 1; i++) {
    ia_[i] += 1;
  }
  for (int i = 0; i < nnz; i++) {
    ja_[i] += 1;
  }

  /* -------------------------------------------------------------------- */
  /*  .. pardiso_chk_matrix(...)                                          */
  /*     Checks the consistency of the given matrix.                      */
  /*     Use this functionality only for debugging purposes               */
  /* -------------------------------------------------------------------- */

  //  pardiso_chkmatrix(&mtype, &n, a, ia, ja, &error);
  //  if (error != 0) {
  //    printf("\nERROR in consistency of matrix: %d", error);
  //    exit(1);
  //  }

  /* -------------------------------------------------------------------- */
  /* ..  pardiso_chkvec(...)                                              */
  /*     Checks the given vectors for infinite and NaN values             */
  /*     Input parameters (see PARDISO user manual for a description):    */
  /*     Use this functionality only for debugging purposes               */
  /* -------------------------------------------------------------------- */

  //  pardiso_chkvec(&n, &nrhs, vec_b, &error);
  //  if (error != 0) {
  //    printf("\nERROR  in right hand side: %d", error);
  //    exit(1);
  //  }

  /* -------------------------------------------------------------------- */
  /* .. pardiso_printstats(...)                                           */
  /*    prints information on the matrix to STDOUT.                       */
  /*    Use this functionality only for debugging purposes                */
  /* -------------------------------------------------------------------- */

  //    pardiso_printstats (&mtype, &n, a, ia, ja, &nrhs, b, &error);
  //  if (error != 0) {
  //    printf("\nERROR right hand side: %d", error);
  //    exit(1);
  //  }

  /* -------------------------------------------------------------------- */
  /* ..  Reordering and Symbolic Factorization.  This step also allocates */
  /*     all memory that is necessary for the factorization.              */
  /* -------------------------------------------------------------------- */
  phase_ = 11;

  pardiso(pt_, &maxfct_, &mnum_, &mtype_, &phase_, &n, a_, ia_, ja_, &idum_,
          &nrhs_, iparm_, &msglvl_, &ddum_, &ddum_, &error_, dparm_);

  if (error_ != 0) {
    printf("\nERROR during symbolic factorization: %d", error_);
    exit(1);
  }
  //    printf("\nReordering completed ... ");
  //    printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
  //    printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

  /* -------------------------------------------------------------------- */
  /* ..  Numerical factorization.                                         */
  /* -------------------------------------------------------------------- */
  phase_ = 22;

  pardiso(pt_, &maxfct_, &mnum_, &mtype_, &phase_, &n, a_, ia_, ja_, &idum_,
          &nrhs_, iparm_, &msglvl_, &ddum_, &ddum_, &error_, dparm_);

  if (error_ != 0) {
    printf("\nERROR during numerical factorization: %d", error_);
    exit(2);
  }
  //    printf("\nFactorization completed ...\n ");

  /* -------------------------------------------------------------------- */
  /* ..  Back substitution and iterative refinement.                      */
  /* -------------------------------------------------------------------- */
  phase_ = 33;

  pardiso(pt_, &maxfct_, &mnum_, &mtype_, &phase_, &n, a_, ia_, ja_, &idum_,
          &nrhs_, iparm_, &msglvl_, vec_b, x_diff, &error_, dparm_);

  if (error_ != 0) {
    printf("\nERROR during solution: %d", error_);
    exit(3);
  }

  //    printf("\nSolve completed ... ");
  //    printf("\nThe solution of the system is: ");
  //    for (i = 0; i < n; i++) {
  //        printf("\n x [%d] = % f", i, x[i] );
  //    }
  //    printf ("\n");
  /* -------------------------------------------------------------------- */
  /* ..  Convert matrix back to 0-based C-notation.                       */
  /* -------------------------------------------------------------------- */
  for (int i = 0; i < n + 1; i++) {
    ia_[i] -= 1;
  }
  for (int i = 0; i < nnz; i++) {
    ja_[i] -= 1;
  }

  /* -------------------------------------------------------------------- */
  /* ..  Termination and release of memory.                               */
  /* -------------------------------------------------------------------- */
  phase_ = -1; /* Release internal memory. */

  pardiso(pt_, &maxfct_, &mnum_, &mtype_, &phase_, &n, &ddum_, ia_, ja_, &idum_,
          &nrhs_, iparm_, &msglvl_, &ddum_, &ddum_, &error_, dparm_);

  // return x_diff
  Eigen::VectorXd x_diff_vec(n);
  for (int i = 0; i < n; i++) {
    x_diff_vec[i] = x_diff[i];
  }

  return x_diff_vec;
}