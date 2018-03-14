//
//  Shell.cpp
//  stap++
//
//  Created by six on 2017/12/9.
//

#include "Shell.h"
//the storage location of the matrixs above is in this header file 
//WE DON'T NEED to consider it anymore in the calculation below
//just use it like in matlab (^_^)
#include "Shell_Def.h" 
//a black magic powered by six!
#include <iostream>
#include <iomanip>
#include <cmath>

//using namespace std;
//template<class T_>
//void debugout(T_* a,int n){
//    for(int i=0;i<n;i++){
//        cout<< a[i]<<endl;
//    }
//}
//
//void debugoutstiff(double*a){
//    for(int i=0;i<12;i++){
//        for(int j=0;j<12;j++){
//            cout<< setiosflags(ios::scientific) <<setprecision(3) << setw(12)<<((i>j)?a[i*(i+3)/2-j]:a[j*(j+3)/2-i]);
//        }
//        cout<< endl;
//    }
//}

// Constructor
CShell::CShell(){
    NEN=9;
    nodes= new CNode*[NEN];

    ND=45;
    LocationMatrix= new unsigned int[ND];

    BmatTotal=new double[5*45*18];
    JacobiDet=new double[18];
    Cttmat=new double[4];

    ElementMaterial = NULL;
}

CShell::~CShell(){
    delete [] nodes;
    delete [] LocationMatrix;
    delete [] BmatTotal;
    delete [] JacobiDet;
    delete [] Cttmat;
}

// Read element data from stream Input
bool CShell::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList){
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
    ElementMaterial = &(dynamic_cast<CShellMaterial*>(MaterialSets))[N - 1];
    
    return true;
}

//	Write element data to stream OutputFile
void CShell::Write(COutputter& output, unsigned int Ele)
{
    output << setw(5) << Ele+1;
    for(int i=0;i<NEN;i++){cout<< setw(9)<< nodes[i]->NodeNumber;}
    cout << setw(12) << ElementMaterial->nset << endl;

}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !

void CShell::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN; N++)
        for (unsigned int D = 0; D < 5; D++)
            LocationMatrix[i++] = nodes[N]->bcode[D];
}

inline void CShell::CalCttmat(){
    double E=(ElementMaterial)->E;
    double posi_rate=(dynamic_cast<CShellMaterial*>(ElementMaterial))->posi_rate;
    Cttmat(1,1)=E/(1-posi_rate*posi_rate);
    Cttmat(1,2)=posi_rate*Cttmat(1,1);
    Cttmat(3,3)=(1-posi_rate)*Cttmat(1,1)/2;
    Cttmat(4,4)=Cttmat(3,3)/1.2;
    Cttmat(5,5)=Cttmat(4,4);
}

void CShell::CalBmatSingle(double* Bmat,double* JacobiDet,double** coor,double yita,double psi,double zeta){

    double h=((CShellMaterial*)ElementMaterial)->thickness/2;
    double N[9];// shape function
    N[0]=psi*(1-psi)/2;             N[1]=1-psi*psi;             N[2]=psi*(1+psi)/2;
    N[3]=N[0]*(1-yita*yita);        N[4]=N[1]*(1-yita*yita);    N[5]=N[2]*(1-yita*yita);
    N[6]=N[0]*yita*(1+yita)/2;      N[7]=N[1]*yita*(1+yita)/2;  N[8]=N[2]*yita*(1+yita)/2;
    N[0]=N[0]*yita*(1-yita)/2;      N[0]=N[1]*yita*(1-yita)/2;  N[0]=N[2]*yita*(1-yita)/2;

    double GN[2][9];//the (psi,yita) Partial derivative of shape function 
    GN(1,1)=0.5-psi;                    GN(1,2)=-2*psi;                     GN(1,3)=0.5+psi;  
    GN(1,4)=-(1-yita*yita)*GN(1,1);     GN(1,5)=-(1-yita*yita)*GN(1,2);     GN(1,6)=-(1-yita*yita)*GN(1,3);
    GN(1,7)=0.5*yita*(1+yita)*GN(1,1);  GN(1,8)=0.5*yita*(1+yita)*GN(1,2);  GN(1,9)=0.5*yita*(1+yita)*GN(1,3);
    GN(1,1)=0.5*yita*(1-yita)*GN(1,1);  GN(1,2)=0.5*yita*(1-yita)*GN(1,2);  GN(1,3)=0.5*yita*(1-yita)*GN(1,3);
    
    GN(2,1)=0.5-yita;                   GN(2,2)=-2*yita;                    GN(2,3)=0.5+yita;  
    GN(2,4)=-(1-psi*psi)*GN(1,1);       GN(2,5)=-(1-psi*psi)*GN(1,2);       GN(2,6)=-(1-psi*psi)*GN(1,3);
    GN(2,7)=0.5*psi*(1+psi)*GN(1,1);    GN(2,8)=0.5*psi*(1+psi)*GN(1,2);    GN(2,9)=0.5*psi*(1+psi)*GN(1,3);
    GN(2,1)=0.5*psi*(1-psi)*GN(1,1);    GN(2,2)=0.5*psi*(1-psi)*GN(1,2);    GN(2,3)=0.5*psi*(1-psi)*GN(1,3);
    

    // considering that J(3,1)=J(3,2)=0 in our situation，the Jacobi_inv canbe easy calculated.
    double Jacobi[3][3];
    for(int j=0;j<3;j++){
        for(int i=0;i<9;i++){
            Jacobi[0][j]+=GN[0][i]*coor[i][j];
            Jacobi[1][j]+=GN[1][i]*coor[i][j];
        }
    }

    (*JacobiDet)=Jacobi(1,1)*Jacobi(2,2)-Jacobi(2,1)*Jacobi(1,2);
    double Jacobi_inv[3][3];
    Jacobi_inv(1,1) =Jacobi(2,2)/(*JacobiDet);   Jacobi_inv(1,2) =-Jacobi(1,2)/(*JacobiDet);
    Jacobi_inv(2,1) =-Jacobi(2,1)/(*JacobiDet);  Jacobi_inv(2,2) =Jacobi(1,1)/(*JacobiDet);
    Jacobi_inv(3,3)=1/h;
    Jacobi_inv(1,3)=-Jacobi_inv(3,3)*(Jacobi_inv(1,1)*Jacobi(1,3)+Jacobi_inv(1,2)*Jacobi(2,3));
    Jacobi_inv(2,3)=-Jacobi_inv(3,3)*(Jacobi_inv(2,1)*Jacobi(1,3)+Jacobi_inv(2,2)*Jacobi(2,3));
    (*JacobiDet)*=h;

    //shapefunction is H*J^(-1) *  [GN(2*9);N(1*9)]**[A(3*5);A(3*5);B(3*5)]
    // ** means convolution
    double A[3][5];
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            A[i][j]=Jacobi_inv[i][j];
        }
    }
    for(int i=0;i<3;i++){A[i][3]=-zeta*Jacobi_inv[i][1]/2;}
    for(int i=0;i<3;i++){A[i][4]=zeta*Jacobi_inv[i][0]/2;}

    double B[3][5];
    for(int i=0;i<3;i++){B[i][3]=-Jacobi_inv[i][1]/2;}
    for(int i=0;i<3;i++){B[i][4]=Jacobi_inv[i][0]/2;}


    int index;
    for(int i=0;i<9;i++){
        // the five rows of Bmat is:
        // Bmat(i*5:45:end)=
        //     GN[0][i]*A(0,:)
        //     GN[1][i]*A(1,:)
        //     GN[0][i]*A(1,:)+GN[1][i]*A(0,:)
        //     GN[1][i]*A(2,:)+N[i]*B(1,:)
        //     GN[0][i]*A(2,:)+N[i]*B(0,:)
        index=5*i;
        for(int item=0;item<5;item++){Bmat[index+item]=GN[0][i]*A[0][item];}
        index=5*i+45;
        for(int item=0;item<5;item++){Bmat[index+item]=GN[1][i]*A[1][item];}
        index=5*i+45*2;
        for(int item=0;item<5;item++){Bmat[index+item]=GN[0][i]*A[1][item]+GN[1][i]*A[0][item];}
        index=5*i+45*3;
        for(int item=0;item<5;item++){Bmat[index+item]=GN[1][i]*A[2][item]+N[i]*B[1][item];}
        index=5*i+45*4;
        for(int item=0;item<5;item++){Bmat[index+item]=GN[0][i]*A[2][item]+N[i]*B[0][item];}
    }

}

inline void CShell::CalBmatTotal(){
    #define Bmat_size 225
    double yita[9]={-0.774596669241483,-0.774596669241483,-0.774596669241483,0,0,0,0.774596669241483,0.774596669241483,0.774596669241483};
    double psi[9]={-0.774596669241483,0,0.774596669241483,-0.774596669241483,0,0.774596669241483,-0.774596669241483,0,0.774596669241483};

    double** coor=new double*[9];
    for(int i=0;i<9;i++){coor[i]=new double[5];}

    for(int i=0;i<9;i++){
        for(int j=0;j<5;j++){
            coor[i][j]=nodes[i]->XYZ[j];
        }
    }

    for(int i=0;i<NEN;i++){CalBmatSingle(BmatTotal+i*Bmat_size,JacobiDet+i,coor,yita[i],psi[i],-1/sqrt(3));}
    
}

inline unsigned int CShell::SizeOfStiffnessMatrix() { return 1035; }

void CShell::ElementStiffnessSingle(double* StiffMatrix,double* Bmat,double* Cttmat,double JacobiDet, double weight){

for(int i=0;i<45;i++){
    for(int j=i;j<45;j++){
        StiffMatrix[j*(j+3)/2-i]+=weight*(
            Bmat[i]*(Cttmat[0]*Bmat[j]+Cttmat[1]*Bmat[j+45])+Bmat[i+45]*(Cttmat[1]*Bmat[j]+Cttmat[0]*Bmat[j+45])
            +Cttmat[2]*Bmat[i+90]*Bmat[j+90]
            +Cttmat[3]*Bmat[i+135]*Bmat[j+135]
            +Cttmat[3]*Bmat[i+180]*Bmat[j+180]
                                          );
    }
}

}

void CShell::ElementStiffness(double* StiffMatrix){
    clear(StiffMatrix, SizeOfStiffnessMatrix());
    CalCttmat();
    CalBmatTotal();
    double weight[9]={0.308641975308642, 0.493827160493827 ,0.308641975308642, 0.493827160493827, 0.790123456790123, 0.493827160493827, 0.308641975308642, 0.493827160493827, 0.308641975308642};
    for(int i=0;i<9;i++){
        ElementStiffnessSingle(StiffMatrix,BmatTotal+i*Bmat_size,Cttmat,JacobiDet[i],weight[i]);
    }
}

void CShell::ElementGravity(double* bodyforce, double Gravity)
{
}

void CShell::ElementStress(double* stress, double* Displacement)
{
//    CShellMaterial* material = dynamic_cast<CShellMaterial*>(ElementMaterial);    // Pointer to material of the element

    //needed to bu improved::3D-case
    double d[45];
    int item=0;
    for(int i=0;i<45;i++){
        d[i]=LocationMatrix[i]?Displacement[LocationMatrix[i]-1]:0;
    }
	
    double sigma[9][5];//目前输出的仍然是高斯点应力
	for (int i = 0; i < 9; i++)
	{
        double* BBBB=BmatTotal+i*Bmat_size;
        for(int j=0;j<5;j++){
            sigma[i][j]=0;
            for(int k=0;k<45;k++){
                sigma[i][j]+=BBBB[45*j+k]*d[k];
            }
        }

        double item=sigma[i][0];
        sigma[i][0]=sigma[i][0]*Cttmat[0]+sigma[i][1]*Cttmat[1];
        sigma[i][1]=item*Cttmat[1]+sigma[i][1]*Cttmat[0];
        sigma[i][2]=sigma[i][2]*Cttmat[2];
        sigma[i][3]=sigma[i][3]*Cttmat[3];
        sigma[i][4]=sigma[i][4]*Cttmat[3];
	}

    for (int i=0;i<4;i++){
        stress[i]=sqrt(sigma[i][0]*sigma[i][0]+sigma[i][1]*sigma[i][1]+4*sigma[i][2]*sigma[i][2]+4*sigma[i][3]*sigma[i][3]+4*sigma[i][4]*sigma[i][4]);
    }
    
}
