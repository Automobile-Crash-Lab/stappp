//
//  Solver_Sparse.cpp
//  stap++
//
//  Created by six on 2017/12/21.
//
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string>
#include "Solver_Sparse.h"


using namespace std;

/* PARDISO prototype. */
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *,
                             double *, int    *,    int *, int *,   int *, int *,
                             int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                                    double *, int *);



Solver_Sparse::Solver_Sparse(){
    n=0;ia=nullptr;ja=nullptr;fabed=false;
}

Solver_Sparse::Solver_Sparse(int _NEQ,int* _ia,int* _ja,double* _a){
    n=_NEQ; ia=_ia; ja=_ja; a=_a;fabed=false;x=new double[n];
    
}

void Solver_Sparse::BackSubstitution(double *b){
    
    /* -------------------------------------------------------------------- */
    /* ..  Back substitution and iterative refinement.                      */
    /* -------------------------------------------------------------------- */
    phase = 33;
    
    iparm[7] = 1;       /* Max numbers of iterative refinement steps. */
    
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, b, x, &error,  dparm);
    
    if (error != 0) {
        printf("\nERROR during solution: %d", error);
        exit(3);
    }
    
    for(int i=0;i<n;i++){b[i]=x[i];}
    

}


void Solver_Sparse::release(){
    
    /* -------------------------------------------------------------------- */
    /* ..  Convert matrix back to 0-based C-notation.                       */
    /* -------------------------------------------------------------------- */
    for (int i = 0; i < n+1; i++) {--ia[i];}
    for (int i = 0; i < nnz; i++) {--ja[i];}
    
    /* -------------------------------------------------------------------- */
    /* ..  Termination and release of memory.                               */
    /* -------------------------------------------------------------------- */
    phase = -1;                 /* Release internal memory. */
    
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
}


void Solver_Sparse::LDLT(){
    nnz = ia[n];
    mtype = -2;        /* Real symmetric matrix */
    nrhs = 1;          /* Number of right hand sides. */
    error = 0;
    solver = 0; /* use sparse direct solver */
    pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);
    
    /* Numbers of processors, value of OMP_NUM_THREADS */
    char* var = getenv("OMP_NUM_THREADS");
    if(var != NULL)
        sscanf( var, "%d", &num_procs );
    else {
        printf("Set environment OMP_NUM_THREADS to 1");
        exit(1);
    }
    iparm[2]  = num_procs;
    
    maxfct = 1;        /* Maximum number of numerical factorizations.  */
    mnum   = 1;         /* Which factorization to use. */
    
    msglvl = 1;         /* Print statistical information  */
    error  = 0;         /* Initialize error flag */
    
    /* -------------------------------------------------------------------- */
    /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
    /*     notation.                                                        */
    /* -------------------------------------------------------------------- */
    for (int i = 0; i < n+1; i++) {
        ++ia[i];
    }
    for (int i = 0; i < nnz; i++) {
        ++ja[i];
    }
    
    /* -------------------------------------------------------------------- */
    /* ..  Reordering and Symbolic Factorization.  This step also allocates */
    /*     all memory that is necessary for the factorization.              */
    /* -------------------------------------------------------------------- */
    phase = 11;

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error, dparm);

    if (error != 0) {
        printf("\nERROR during symbolic factorization: %d", error);
        exit(1);
    }
    printf("\nReordering completed ... ");
    printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
    printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

    /* -------------------------------------------------------------------- */
    /* ..  Numerical factorization.                                         */
    /* -------------------------------------------------------------------- */
    phase = 22;
    iparm[32] = 1; /* compute determinant */

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);

    if (error != 0) {
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }
    printf("\nFactorization completed ...\n ");
    
}


void Solver_Sparse::solve(double* b){
    nnz = ia[n];
    mtype = 2;        /* Real symmetric matrix */
    nrhs = 1;          /* Number of right hand sides. */
    
    /* -------------------------------------------------------------------- */
    /* ..  Setup Pardiso control parameters.                                */
    /* -------------------------------------------------------------------- */
    
    error = 0;
    solver = 0; /* use sparse direct solver */
    pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);
    
    if (error != 0)
    {
        if (error == -10 )
            printf("No license file found \n");
        if (error == -11 )
            printf("License is expired \n");
        if (error == -12 )
            printf("Wrong username or hostname \n");
    }
    else
        printf("[PARDISO]: License check was successful ... \n");
    
    /* Numbers of processors, value of OMP_NUM_THREADS */
    char* var = getenv("OMP_NUM_THREADS");
    if(var != NULL)
        sscanf( var, "%d", &num_procs );
    
    else {
        printf("Set environment OMP_NUM_THREADS to 1");
        exit(1);
    }
    iparm[2]  = num_procs;
    
    maxfct = 1;        /* Maximum number of numerical factorizations.  */
    mnum   = 1;         /* Which factorization to use. */
    
    msglvl = 1;         /* Print statistical information  */
    error  = 0;         /* Initialize error flag */
    
    /* -------------------------------------------------------------------- */
    /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
    /*     notation.                                                        */
    /* -------------------------------------------------------------------- */
    for (int i = 0; i < n+1; i++) {
        ++ia[i];
    }
    for (int i = 0; i < nnz; i++) {
        ++ja[i];
    }
    
    /* -------------------------------------------------------------------- */
    /*  .. pardiso_chk_matrix(...)                                          */
    /*     Checks the consistency of the given matrix.                      */
    /*     Use this functionality only for debugging purposes               */
    /* -------------------------------------------------------------------- */
    
//    pardiso_chkmatrix  (&mtype, &n, a, ia, ja, &error);
//    if (error != 0) {
//        printf("\nERROR in consistency of matrix: %d", error);
//        exit(1);
//    }
    
    /* -------------------------------------------------------------------- */
    /* ..  pardiso_chkvec(...)                                              */
    /*     Checks the given vectors for infinite and NaN values             */
    /*     Input parameters (see PARDISO user manual for a description):    */
    /*     Use this functionality only for debugging purposes               */
    /* -------------------------------------------------------------------- */
    
//    pardiso_chkvec (&n, &nrhs, b, &error);
//    if (error != 0) {
//        printf("\nERROR  in right hand side: %d", error);
//        exit(1);
//    }
//
    /* -------------------------------------------------------------------- */
    /* .. pardiso_printstats(...)                                           */
    /*    prints information on the matrix to STDOUT.                       */
    /*    Use this functionality only for debugging purposes                */
    /* -------------------------------------------------------------------- */
    
//    pardiso_printstats (&mtype, &n, a, ia, ja, &nrhs, b, &error);
//    if (error != 0) {
//        printf("\nERROR right hand side: %d", error);
//        exit(1);
//    }
    
    /* -------------------------------------------------------------------- */
    /* ..  Reordering and Symbolic Factorization.  This step also allocates */
    /*     all memory that is necessary for the factorization.              */
    /* -------------------------------------------------------------------- */
    phase = 11;
    
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error, dparm);
    
    if (error != 0) {
        printf("\nERROR during symbolic factorization: %d", error);
        exit(1);
    }
    printf("\nReordering completed ... ");
    printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
    printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
    
    /* -------------------------------------------------------------------- */
    /* ..  Numerical factorization.                                         */
    /* -------------------------------------------------------------------- */
    phase = 22;
    iparm[32] = 1; /* compute determinant */
    
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
    
    if (error != 0) {
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }
    printf("\nFactorization completed ...\n ");
    
    /* -------------------------------------------------------------------- */
    /* ..  Back substitution and iterative refinement.                      */
    /* -------------------------------------------------------------------- */
    phase = 33;
    
    iparm[7] = 1;       /* Max numbers of iterative refinement steps. */
    
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, b, x, &error,  dparm);
    
    if (error != 0) {
        printf("\nERROR during solution: %d", error);
        exit(3);
    }
    
    printf("\nSolve completed ... ");
    printf("\nThe solution of the system is: ");
    for (int i = 0; i < n; i++) {
        printf("\n x [%d] = % f", i, x[i] );
    }
    printf ("\n");
    
    /* -------------------------------------------------------------------- */
    /* ..  Convert matrix back to 0-based C-notation.                       */
    /* -------------------------------------------------------------------- */
    for (int i = 0; i < n+1; i++) {
        --ia[i];
    }
    for (int i = 0; i < nnz; i++) {
        --ja[i];
    }
    
    /* -------------------------------------------------------------------- */
    /* ..  Termination and release of memory.                               */
    /* -------------------------------------------------------------------- */
    phase = -1;                 /* Release internal memory. */
    
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
    
}


//#include <mkl.h>
////pardisoinit (_MKL_DSS_HANDLE_t pt, MKL_INT *mtype, MKL_INT *iparm);
//void Solver_Sparse::LDLT(){
//    nnz = ia[n];
//    mtype = -2;        /* Real symmetric matrix */
//    nrhs = 1;          /* Number of right hand sides. */
//    error = 0;
//    solver = 0; /* use sparse direct solver */
//    pardisoinit (pt,  &mtype, iparm);
//
//    maxfct = 1;        /* Maximum number of numerical factorizations.  */
//    mnum   = 1;         /* Which factorization to use. */
//
//    msglvl = 0;         /* Print statistical information  */
//    iparm[35]=1; //ia and ja is zero-based
//
//    /* -------------------------------------------------------------------- */
//    /* ..  Reordering and Symbolic Factorization.  This step also allocates */
//    /*     all memory that is necessary for the factorization.              */
//    /* -------------------------------------------------------------------- */
//    phase = 11;
//
//    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
//             &n, a, ia, ja, &idum, &nrhs,
//             iparm, &msglvl, &ddum, &ddum, &error, dparm);
//
//    if (error != 0) {
//        printf("\nERROR during symbolic factorization: %d", error);
//        exit(1);
//    }
//    printf("\nReordering completed ... ");
//    printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
//    printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
//
//    /* -------------------------------------------------------------------- */
//    /* ..  Numerical factorization.                                         */
//    /* -------------------------------------------------------------------- */
//    phase = 22;
//    iparm[32] = 1; /* compute determinant */
//
//    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
//             &n, a, ia, ja, &idum, &nrhs,
//             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
//
//    if (error != 0) {
//        printf("\nERROR during numerical factorization: %d", error);
//        exit(2);
//    }
//    printf("\nFactorization completed ...\n ");
//
//}


