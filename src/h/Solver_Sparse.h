//
//  Solver_Sparse.h
//  stap++
//
//  Created by six on 2017/12/21.
//

#ifndef Solver_Sparse_h
#define Solver_Sparse_h
//#include "SparseMatrix.h"
class Solver_Sparse{
public:
    int n;
    int*ia,*ja;
    double *a;
    double *x;
    bool fabed;
//pardiso needed
    int      nnz;
    int      mtype;        /* Real symmetric matrix */
    int      nrhs;          /* Number of right hand sides. */
    
    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */
    void    *pt[64];
    
    /* Pardiso control parameters. */
    int      iparm[64];
    double   dparm[64];
    int      maxfct, mnum, phase, error, msglvl, solver;
    
    /* Number of processors. */
    int      num_procs;
    
    /* Auxiliary variables. */

    double   ddum;              /* Double dummy */
    int      idum;              /* Integer dummy. */
//pardiso needed
    
    Solver_Sparse();
    Solver_Sparse(int,int*,int*,double*);
    void LDLT();
    void BackSubstitution(double* Force);
    void release();
    void solve(double* b);
    
};


#endif /* Solver_Sparse_h */
