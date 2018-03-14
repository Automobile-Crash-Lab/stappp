//
//  CSparseMatrix.h
//  stap++
//
//  Created by six on 2017/12/21.
//

#ifndef CSparseMatrix_h
#define CSparseMatrix_h
class CSparseMatrix{
public:
    int NEQ;
    int*ia,*ja;
    double *a;
    Solver_Sparse();
    Solver_Sparse(int,int*,int*,double*);
    void BackSubstitution(double* Force);
}

CSparseMatrix::CSparseMatrix(){
    NEQ=0;ia=nullptr;ja=nullptr;a=nullptr;
}

CSparseMatrix::CSparseMatrix(int _NEQ;int* _ia,int* _ja,double* _a){
    NEQ=_NEQ; ia=_ia; ja=_ja; a=_a;
}

#endif /* CSparseMatrix_h */
