//
//  4Q.cpp
//  stap++
//
//  Created by six on 2017/11/9.
//

#include "4Q.h"
//the storage location of the matrixs above is in this header file 
//WE DON'T NEED to consider it anymore in the calculation below
//just use it like in matlab (^_^)
#include "4Q_Def.h" 
//a black magic powered by six!
#include <iostream>
#include <iomanip>
#include <cmath>

// using namespace std;
// template<class T_>
// void debugout(T_* a,int n){
//     for(int i=0;i<n;i++){
//         cout<< a[i]<<endl;
//     }
// }

// void debugoutstiff(double*a){
//     for(int i=0;i<12;i++){
//         for(int j=0;j<12;j++){
//             cout<< setiosflags(ios::scientific) <<setprecision(3) << setw(12)<<((i>j)?a[i*(i+3)/2-j]:a[j*(j+3)/2-i]);
//         }
//         cout<< endl;
//     }
// }

// Constructor
C4Q::C4Q(){
    NEN=4;
    nodes= new CNode*[NEN];

    ND=12;
    LocationMatrix= new unsigned int[ND];

    ElementMaterial = NULL;
}

C4Q::~C4Q(){
    delete [] nodes;
    delete [] LocationMatrix;
}

// Read element data from stream Input
bool C4Q::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList){
    unsigned int N;

    Input >> N;

    if (N != Ele + 1){
		cout << "*** Error *** Elements must be inputted in order !" << endl 
			 << "   Expected element : " << Ele + 1 << endl
			 << "   Provided element : " << N << endl;
		return false;
    }
    
    // N turns to the number of node
    for(int i=0;i<NEN;i++){
        Input >> N;
        nodes[i]=&NodeList[N - 1];
    }

    //N turns to material property set number
    Input >> N;
    ElementMaterial = &(dynamic_cast<C4QMaterial*>(MaterialSets))[N - 1];
    
    return true;
}

//	Write element data to stream OutputFile
void C4Q::Write(COutputter& output, unsigned int Ele)
{
    output << setw(5) << Ele+1;
    for(int i=0;i<NEN;i++){cout<< setw(9)<< nodes[i]->NodeNumber;}
    cout << setw(12) << ElementMaterial->nset << endl;

}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !

void C4Q::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN; N++)
        for (unsigned int D = 0; D < 3; D++)
            LocationMatrix[i++] = nodes[N]->bcode[D];
    
}

inline void C4Q::CalCttmat(){
    double E=(ElementMaterial)->E;
    double posi_rate=(dynamic_cast<C4QMaterial*>(ElementMaterial))->posi_rate;
    Cttmat(1,1)=E/(1-posi_rate*posi_rate);
    Cttmat(1,2)=posi_rate*Cttmat(1,1);
    Cttmat(3,3)=(1-posi_rate)*Cttmat(1,1)/2;
    
    
    
}

void C4Q::CalBmatSingle(double* Bmat,double* JacobiDet,double* coor,double yita,double psi){
    double GN[8];
    GN(1,1)=0.25*(yita-1);  GN(1,2)=-GN(1,1);       GN(1,3)=0.25*(yita+1);  GN(1,4)=-GN(1,3);
    GN(2,1)=0.25*(psi-1);   GN(2,2)=-0.25*(psi+1);  GN(2,3)=-GN(2,2);       GN(2,4)=-GN(2,1);
    
    double Jacobi[4];
#define CalJacobi(i,j) (Jacobi(i,j)=GN(i,1) *coor(1,j) + GN(i,2) *coor(2,j) + GN(i,3) *coor(3,j) + GN(i,4) *coor(4,j))
    CalJacobi(1,1); CalJacobi(1,2);
    CalJacobi(2,1); CalJacobi(2,2);
    
    (*JacobiDet)=Jacobi(1,1)*Jacobi(2,2)-Jacobi(2,1)*Jacobi(1,2);
    
    double Jacobi_inv[4];
    Jacobi_inv(1,1) =Jacobi(2,2)/(*JacobiDet);   Jacobi_inv(1,2) =-Jacobi(1,2)/(*JacobiDet);
    Jacobi_inv(2,1) =-Jacobi(2,1)/(*JacobiDet);  Jacobi_inv(2,2) =Jacobi(1,1)/(*JacobiDet);
    
#define CalBmat(i,j) (Bmat(i,j)=Jacobi_inv(i,1) *GN(1,j) + Jacobi_inv(i,2) *GN(2,j))
    CalBmat(1,1);   CalBmat(1,2);   CalBmat(1,3);   CalBmat(1,4);
    CalBmat(2,1);   CalBmat(2,2);   CalBmat(2,3);   CalBmat(2,4);
}

inline void C4Q::CalBmatTotal(){
    #define Bmat_size 8
    double item=1/sqrt(3);
    double yita[4]={-item,-item,item,item};
    double psi[4]={-item,item,-item,item};
    double coor[8];
    coor(1,1)=nodes[0]->XYZ[0];  coor(1,2)=nodes[0]->XYZ[1];
    coor(2,1)=nodes[1]->XYZ[0];  coor(2,2)=nodes[1]->XYZ[1];
    coor(3,1)=nodes[2]->XYZ[0];  coor(3,2)=nodes[2]->XYZ[1];
    coor(4,1)=nodes[3]->XYZ[0];  coor(4,2)=nodes[3]->XYZ[1];

    for(int i=0;i<4;i++){CalBmatSingle(BmatTotal+i*Bmat_size,JacobiDet+i,coor,yita[i],psi[i]);}

}

inline unsigned int C4Q::SizeOfStiffnessMatrix() { return 78; }

void C4Q::ElementStiffnessSingle(double* StiffMatrix,double* Bmat,double* Cttmat,double JacobiDet, double weight){
    CalCttmat();
#define CalStiff(i,j) (StiffMatrix(i,j)+=weight*JacobiDet*(NBmat(1,i)*(Cttmat(1,1)*NBmat(1,j)+Cttmat(1,2)*NBmat(2,j))+NBmat(2,i)*(Cttmat(2,1)*NBmat(1,j)+Cttmat(2,2)*NBmat(2,j))+NBmat(3,i)*Cttmat(3,3)*NBmat(3,j)))
    CalStiff(1,1);
    CalStiff(2,1);CalStiff(2,2);
    CalStiff(3,1);CalStiff(3,2);CalStiff(3,3);
    CalStiff(4,1);CalStiff(4,2);CalStiff(4,3);CalStiff(4,4);
    CalStiff(5,1);CalStiff(5,2);CalStiff(5,3);CalStiff(5,4);CalStiff(5,5);
    CalStiff(6,1);CalStiff(6,2);CalStiff(6,3);CalStiff(6,4);CalStiff(6,5);CalStiff(6,6);
    CalStiff(7,1);CalStiff(7,2);CalStiff(7,3);CalStiff(7,4);CalStiff(7,5);CalStiff(7,6);CalStiff(7,7);
    CalStiff(8,1);CalStiff(8,2);CalStiff(8,3);CalStiff(8,4);CalStiff(8,5);CalStiff(8,6);CalStiff(8,7);CalStiff(8,8);
}

void C4Q::ElementStiffness(double* StiffMatrix){
    clear(StiffMatrix, SizeOfStiffnessMatrix());
    CalCttmat();
    CalBmatTotal();
    ElementStiffnessSingle(StiffMatrix,BmatTotal,Cttmat,JacobiDet[0],1);
    ElementStiffnessSingle(StiffMatrix,BmatTotal+8,Cttmat,JacobiDet[1],1);
    ElementStiffnessSingle(StiffMatrix,BmatTotal+16,Cttmat,JacobiDet[2],1);
    ElementStiffnessSingle(StiffMatrix,BmatTotal+24,Cttmat,JacobiDet[3],1);
    
}

void C4Q::ElementGravity(double* bodyforce, double Gravity)
{
    clear(bodyforce,12);
}

void C4Q::ElementStress(double* stress, double* Displacement)
{
//    C4QMaterial* material = dynamic_cast<C4QMaterial*>(ElementMaterial);    // Pointer to material of the element

    //needed to bu improved::3D-case
    double d[8];
    int item=0;
    for(int i=0;i<12;i++){
        if(LocationMatrix[i])d[item++]=Displacement[LocationMatrix[i]-1];
    }
	
    double sigma[4][3];
	for (int i = 0; i < 4; i++)
	{
        double* BBBB=BmatTotal+i*Bmat_size;
        sigma[i][0]=BBBB[0]*d[0]+BBBB[1]*d[2]+BBBB[2]*d[4]+BBBB[3]*d[6];
        sigma[i][1]=BBBB[4]*d[1]+BBBB[5]*d[3]+BBBB[6]*d[5]+BBBB[7]*d[7];
        sigma[i][2]=BBBB[4]*d[0]+BBBB[0]*d[1]+BBBB[5]*d[2]+BBBB[1]*d[3]+BBBB[6]*d[4]+BBBB[2]*d[5]+BBBB[7]*d[6]+BBBB[3]*d[7];
        double item=sigma[i][0];
        sigma[i][0]=sigma[i][0]*Cttmat[0]+sigma[i][1]*Cttmat[1];
        sigma[i][1]=item*Cttmat[1]+sigma[i][1]*Cttmat[0];
        sigma[i][2]=sigma[i][2]*Cttmat[2];
	}

    double a1=1+sqrt(3)/2,a2=1-sqrt(3)/2,b=-0.5;
//    a1,b,a2,b
//    b,a1,b,a2
//    a2,b,a1,b
//    b,a2,b,a1
    double sigma2[4][3];
    for (int j=0;j<3;j++){
        sigma2[0][j]=a1*sigma[0][j]+b*sigma[1][j]+a2*sigma[2][j]+b*sigma[3][j];
        sigma2[1][j]=b*sigma[0][j]+a1*sigma[1][j]+b*sigma[2][j]+a2*sigma[3][j];
        sigma2[2][j]=a2*sigma[0][j]+b*sigma[1][j]+a1*sigma[2][j]+b*sigma[3][j];
        sigma2[3][j]=b*sigma[0][j]+a2*sigma[1][j]+b*sigma[2][j]+a1*sigma[3][j];
    }

    for (int i=0;i<4;i++){
        stress[i]=sqrt(sigma2[i][0]*sigma2[i][0]+sigma2[i][1]*sigma2[i][1]+4*sigma2[i][2]*sigma2[i][2]);
    }
    
}
