#include "pardiso_unsym.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

/* PARDISO prototype. */
extern "C" void pardisoinit(void*, int*, int*, int*, double*, int*);
extern "C" void pardiso(void*, int*, int*, int*, int*, int*, double*, int*,
                        int*, int*, int*, int*, int*, double*, double*, int*,
                        double*);
extern "C" void pardiso_chkmatrix(int*, int*, double*, int*, int*, int*);
extern "C" void pardiso_chkvec(int*, int*, double*, int*);
extern "C" void pardiso_printstats(int*, int*, double*, int*, int*, int*,
                                   double*, int*);

bool pipenetwork::Pardiso_unsym::solve() {
  // configure pardiso
  int n = this->vec_b_->size();
  double* vec_x = this->vec_x_->data();
  double* vec_b = this->vec_b_->data();
  double x_diff[n];

  auto ia = (this->mat_a_)->outerIndexPtr();
  auto ja = (this->mat_a_)->innerIndexPtr();
  auto a = (this->mat_a_)->valuePtr();

  int nnz = ia[n];
  int mtype = 11; /* Real unsymmetric matrix */

  /* RHS and solution vectors. */
  int nrhs = 1; /* Number of right hand sides. */

  /* Internal solver memory pointer pt,                  */
  /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
  /* or void *pt[64] should be OK on both architectures  */
  void* pt[64];

  /* Pardiso control parameters. */
  int iparm[64];
  double dparm[64];
  int solver;
  int maxfct, mnum, phase, error, msglvl;

  /* Number of processors. */
  int num_procs;

  /* Auxiliary variables. */
  char* var;
  int i;

  double ddum; /* Double dummy */
  int idum;    /* Integer dummy. */

  /* -------------------------------------------------------------------- */
  /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
  /*     notation.                                                        */
  /* -------------------------------------------------------------------- */
  for (i = 0; i < n + 1; i++) {
    ia[i] += 1;
  }
  for (i = 0; i < nnz; i++) {
    ja[i] += 1;
  }

  /* -------------------------------------------------------------------- */
  /* ..  Setup Pardiso control parameters and initialize the solvers      */
  /*     internal adress pointers. This is only necessary for the FIRST   */
  /*     call of the PARDISO solver.                                      */
  /* ---------------------------------------------------------------------*/

  error = 0;
  solver = 0; /* use sparse direct solver */
  pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);

  if (error != 0) {
    if (error == -10) printf("No license file found \n");
    if (error == -11) printf("License is expired \n");
    if (error == -12) printf("Wrong username or hostname \n");
    return 1;
  }
  //    else
  //        printf("[PARDISO]: License check was successful ... \n");

  /* Numbers of processors, value of OMP_NUM_THREADS */
  var = getenv("OMP_NUM_THREADS");
  if (var != NULL)
    sscanf(var, "%d", &num_procs);
  else {
    printf("Set environment OMP_NUM_THREADS to 1");
    exit(1);
  }

  /* -------------------------------------------------------------------- */
  /* .. Setup Pardiso control parameters. */
  /* -------------------------------------------------------------------- */
  for (i = 0; i < 64; i++) {
    iparm[i] = 0;
  }
  iparm[0] = 1; /* No solver default */
  iparm[1] = 0; /* Fill-in reordering from METIS */
  /* Numbers of processors, value of OMP_NUM_THREADS */
  iparm[2] = num_procs;
  iparm[3] = 0;   /* No iterative-direct algorithm */
  iparm[4] = 0;   /* No user fill-in reducing permutation */
  iparm[5] = 0;   /* Write solution into x */
  iparm[6] = 0;   /* Not in use */
  iparm[7] = 2;   /* Max numbers of iterative refinement steps */
  iparm[8] = 0;   /* Not in use */
  iparm[9] = 10;  /* Perturb the pivot elements with 1E-13 */
  iparm[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
  iparm[11] = 0;  /* Not in use */
  iparm[12] = 0;  /* Not in use */
  iparm[13] = 0;  /* Output: Number of perturbed pivots */
  iparm[14] = 0;  /* Not in use */
  iparm[15] = 0;  /* Not in use */
  iparm[16] = 0;  /* Not in use */
  iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
  iparm[18] = -1; /* Output: Mflops for LU factorization */
  iparm[19] = 0;  /* Output: Numbers of CG Iterations */
  maxfct = 1;     /* Maximum number of numerical factorizations. */
  mnum = 1;       /* Which factorization to use. */
  msglvl = 1;     /* Print statistical information in file */
  error = 0;      /* Initialize error flag */

  maxfct = 1; /* Maximum number of numerical factorizations.  */
  mnum = 1;   /* Which factorization to use. */

  msglvl = 0; /* Print statistical information  */
  error = 0;  /* Initialize error flag */

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
  phase = 11;

  pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
          iparm, &msglvl, &ddum, &ddum, &error, dparm);

  if (error != 0) {
    printf("\nERROR during symbolic factorization: %d", error);
    exit(1);
  }
  //    printf("\nReordering completed ... ");
  //    printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
  //    printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

  /* -------------------------------------------------------------------- */
  /* ..  Numerical factorization.                                         */
  /* -------------------------------------------------------------------- */
  phase = 22;

  pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
          iparm, &msglvl, &ddum, &ddum, &error, dparm);

  if (error != 0) {
    printf("\nERROR during numerical factorization: %d", error);
    exit(2);
  }
  //    printf("\nFactorization completed ...\n ");

  /* -------------------------------------------------------------------- */
  /* ..  Back substitution and iterative refinement.                      */
  /* -------------------------------------------------------------------- */
  phase = 33;

  pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
          iparm, &msglvl, vec_b, x_diff, &error, dparm);

  if (error != 0) {
    printf("\nERROR during solution: %d", error);
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
  for (i = 0; i < n + 1; i++) {
    ia[i] -= 1;
  }
  for (i = 0; i < nnz; i++) {
    ja[i] -= 1;
  }

  /* -------------------------------------------------------------------- */
  /* ..  Termination and release of memory.                               */
  /* -------------------------------------------------------------------- */
  phase = -1; /* Release internal memory. */

  pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs,
          iparm, &msglvl, &ddum, &ddum, &error, dparm);

  for (i = 0; i < n; i++) {
    (*vec_x_)[i] -= x_diff[i];
  }

  return true;
}