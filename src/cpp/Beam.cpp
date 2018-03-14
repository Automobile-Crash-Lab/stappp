/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.0, October 14, 2017                                         */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Beam.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdlib.h> 
using namespace std;

//	Constructor
CBeam::CBeam()
{

	NEN = 2;	// Each element has 2 nodes
	nodes = new CNode*[NEN];
    
    ND = 12;
    LocationMatrix = new unsigned int[ND];

	ElementMaterial = NULL;
}

//	Desconstructor
CBeam::~CBeam()
{
	delete [] nodes;
    delete [] LocationMatrix;
}

//	Read element data from stream Input
bool CBeam::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int N;
	
	Input >> N;	// element number

	if (N != Ele + 1)
	{
		cout << "*** Error *** Elements must be inputted in order !" << endl 
			 << "   Expected element : " << Ele + 1 << endl
			 << "   Provided element : " << N << endl;

		return false;
	}

	unsigned int MSet;	// Material property set number
	unsigned int N1, N2; 	// Left node number and right node number

	Input >> N1 >> N2 >>MSet ;

	ElementMaterial = &((CBeamMaterial*)MaterialSets)[MSet - 1];
	nodes[0] = &NodeList[N1 - 1];
	nodes[1] = &NodeList[N2 - 1];
	

	return true;
}

//	Write element data to stream OutputFile
void CBeam::Write(COutputter& output , unsigned int Ele)
{
	output << setw(5) << Ele+1 << setw(11) << nodes[0]->NodeNumber
		 << setw(9) << nodes[1]->NodeNumber << setw(12) << ElementMaterial->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CBeam::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN; N++)
        for (unsigned int D = 0; D < 6; D++)
            LocationMatrix[i++] = nodes[N]->bcode[D];
}


//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 2 node Beam element, element stiffness is a 12x12 matrix, whose upper triangular part
//	has 78 elements
unsigned int CBeam::SizeOfStiffnessMatrix() { return 78; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CBeam::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());
	
//	Calculate element stiffness matrix
	
	
	CBeamMaterial* material = (CBeamMaterial*)ElementMaterial;	// Pointer to material of the element

	//��Ԫ�ն���δת������ϵ��
    double Ke[12][12]={0};
	double k=material-> E*material->Area;
	double le=sqrt((nodes[0]->XYZ[0]-nodes[1]->XYZ[0])*(nodes[0]->XYZ[0]-nodes[1]->XYZ[0])
		      +(nodes[0]->XYZ[1]-nodes[1]->XYZ[1])*(nodes[0]->XYZ[1]-nodes[1]->XYZ[1])
			  +(nodes[0]->XYZ[2]-nodes[1]->XYZ[2])*(nodes[0]->XYZ[2]-nodes[1]->XYZ[2]));
	Ke[0][0]=k/le;
	Ke[0][6]=Ke[6][0]=-k/le;

	Ke[1][1]=12*material-> E*material-> Iz/(le*le*le);
	Ke[1][5]=Ke[5][1]=6*material-> E*material-> Iz/(le*le);
	Ke[1][7]=Ke[7][1]=-12*material-> E*material-> Iz/(le*le*le);
	Ke[1][11]=Ke[11][1]=6*material-> E*material-> Iz/(le*le);

	Ke[2][2]=12*material-> E*material-> Iy/(le*le*le);
	Ke[2][4]=Ke[4][2]=-6*material-> E*material-> Iy/(le*le);
	Ke[2][8]=Ke[8][2]=-12*material-> E*material-> Iy/(le*le*le);
	Ke[2][10]=Ke[10][2]=6*material-> E*material-> Iy/(le*le);

	Ke[3][3]=material-> E*material->Ip/(2*(1+material->v)*le);
	Ke[3][9]=Ke[9][3]=-material-> E*material->Ip/(2*(1+material->v)*le);

	Ke[4][4]=4*material-> E*material-> Iy/(le);
	Ke[4][8]=Ke[8][4]=-6*material-> E*material-> Iy/(le*le);
	Ke[4][10]=Ke[10][4]=2*material-> E*material-> Iy/(le);

	Ke[5][5]=4*material-> E*material-> Iz/(le);
	Ke[5][7]=Ke[7][5]=-6*material-> E*material-> Iz/(le*le);
	Ke[5][11]=Ke[11][5]=2*material-> E*material-> Iz/(le);

	Ke[6][6]=k/le;
	
	Ke[7][7]=12*material-> E*material-> Iz/(le*le*le);
	Ke[7][11]=Ke[11][7]=-6*material-> E*material-> Iz/(le*le);

	Ke[8][8]=12*material-> E*material-> Iy/(le*le*le);
	Ke[8][10]=Ke[10][8]=-6*material-> E*material-> Iy/(le*le);

	Ke[9][9]=material-> E*material->Ip/(2*(1+material->v)*le);

	Ke[10][10]=4*material-> E*material-> Iy/(le);

	Ke[11][11]=4*material-> E*material-> Iz/(le);

    double L[12][12]={0};
	L[0][0]=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0])/le;
	L[0][2]=(nodes[1]->XYZ[2]-nodes[0]->XYZ[2])/le;
	L[2][0]=-(nodes[1]->XYZ[2]-nodes[0]->XYZ[2])/le;
	L[2][2]=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0])/le;
	L[1][1]=1;
	L[3][3]=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0])/le;
	L[3][5]=(nodes[1]->XYZ[2]-nodes[0]->XYZ[2])/le;
	L[5][3]=-(nodes[1]->XYZ[2]-nodes[0]->XYZ[2])/le;
	L[5][5]=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0])/le;
	L[4][4]=1;
	L[6][6]=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0])/le;
	L[6][8]=(nodes[1]->XYZ[2]-nodes[0]->XYZ[2])/le;
	L[8][6]=-(nodes[1]->XYZ[2]-nodes[0]->XYZ[2])/le;
	L[8][8]=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0])/le;
	L[7][7]=1;
	L[9][9]=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0])/le;
	L[9][11]=(nodes[1]->XYZ[2]-nodes[0]->XYZ[2])/le;
	L[11][9]=-(nodes[1]->XYZ[2]-nodes[0]->XYZ[2])/le;
	L[11][11]=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0])/le;
	L[10][10]=1;
	
	
	//ת����ȫ������ϵ
   
	double K[12][12]={0};

	for(unsigned int i = 0; i < 12; i++)
	{
		//cout << endl;
		for(unsigned int j = 0; j < 12; j++)
		{
			for(unsigned int p = 0; p < 12; p++)
	        {
		         for(unsigned int q = 0; q < 12; q++)
		         {
			         K[i][j]+=L[q][i]*Ke[q][p]*L[p][j];
					 
				 }
			}
			//cout << K[i][j]-Ke[i][j];
		}
	}

	
		
		for (unsigned int j = 0; j < 12; j++)
		{
			for (unsigned int i = 0; i <j+1; i++)
			{				
					
				Matrix[(j+1)*(j+2)/2-i-1]=K[i][j];
									
			}
		}
	
}

void CBeam::ElementGravity(double* bodyforce, double Gravity)
{
	clear(bodyforce,12);

// Elemenet material
	CBeamMaterial* material_ = (CBeamMaterial*)ElementMaterial;	// Pointer to material of the element

	double x[2],y[2],z[2]; // Coordinates of 2 nodes
	for (unsigned int i = 0; i < 2; i++)
		{
			x[i]=nodes[i]->XYZ[0];
			y[i]=nodes[i]->XYZ[1];
			z[i]=nodes[i]->XYZ[2];
		}
	double l=sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0])+(z[1]-z[0])*(z[1]-z[0]));
	
	double cos=(x[1]-x[0])/l;
	double q = material_->rho * Gravity;
	bodyforce[0] = 0;
	bodyforce[1] = 0;
	bodyforce[2] = q*l/2;
	bodyforce[3] = 0;
	bodyforce[4] = q*l*l/12;
	bodyforce[5] = 0;
	bodyforce[6] = 0;
	bodyforce[7] = 0;
	bodyforce[8] = q*l/2;
	bodyforce[9] = 0;
	bodyforce[10] = -q*l*l/12;
	bodyforce[11] = 0;
}

//	Calculate element stress 
void CBeam::ElementStress(double* stress, double* Displacement)
{
	CBeamMaterial* material = (CBeamMaterial*)ElementMaterial;	// Pointer to material of the element
	
	double DX[3];	//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	double L2 = 0;	//	Square of bar length (L^2)

	for (unsigned int i = 0; i < 3; i++)
	{
		DX[i] = nodes[1]->XYZ[i] - nodes[0]->XYZ[i];
		L2 = L2 + DX[i]*DX[i];
	}

	double S[6];
	for (unsigned int i = 0; i < 3; i++)
	{
		S[i] = -DX[i] * material->E / L2;
		S[i+3] = -S[i];
	}
	
	stress[0] = 0.0;
	for (unsigned int i = 0; i < 3; i++)
	{
		if (LocationMatrix[i])
			stress[0] += S[i] * Displacement[LocationMatrix[i]-1];
	}
	for (unsigned int i = 6; i < 9; i++)
	{
		if (LocationMatrix[i])
			stress[0] += S[i-3] * Displacement[LocationMatrix[i]-1];
	}
	double d[12]={0};

	double le=sqrt((nodes[0]->XYZ[0]-nodes[1]->XYZ[0])*(nodes[0]->XYZ[0]-nodes[1]->XYZ[0])
		      +(nodes[0]->XYZ[1]-nodes[1]->XYZ[1])*(nodes[0]->XYZ[1]-nodes[1]->XYZ[1])
			  +(nodes[0]->XYZ[2]-nodes[1]->XYZ[2])*(nodes[0]->XYZ[2]-nodes[1]->XYZ[2]));

	double L[12][12]={0};
	L[0][0]=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0])/le;
	L[0][2]=(nodes[1]->XYZ[2]-nodes[0]->XYZ[2])/le;
	L[2][0]=-(nodes[1]->XYZ[2]-nodes[0]->XYZ[2])/le;
	L[2][2]=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0])/le;
	L[1][1]=1;
	L[3][3]=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0])/le;
	L[3][5]=(nodes[1]->XYZ[2]-nodes[0]->XYZ[2])/le;
	L[5][3]=-(nodes[1]->XYZ[2]-nodes[0]->XYZ[2])/le;
	L[5][5]=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0])/le;
	L[4][4]=1;
	L[6][6]=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0])/le;
	L[6][8]=(nodes[1]->XYZ[2]-nodes[0]->XYZ[2])/le;
	L[8][6]=-(nodes[1]->XYZ[2]-nodes[0]->XYZ[2])/le;
	L[8][8]=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0])/le;
	L[7][7]=1;
	L[9][9]=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0])/le;
	L[9][11]=(nodes[1]->XYZ[2]-nodes[0]->XYZ[2])/le;
	L[11][9]=-(nodes[1]->XYZ[2]-nodes[0]->XYZ[2])/le;
	L[11][11]=(nodes[1]->XYZ[0]-nodes[0]->XYZ[0])/le;
	L[10][10]=1;

	for (unsigned int i = 0; i < 12; i++)
		{
			for(unsigned int j = 0; j < 12; j++)
			{
				if(LocationMatrix[j]>0)
				{
					d[i]+=L[i][j]*Displacement[LocationMatrix[j]-1];
				}
			}
		}
	
	stress[1]=material-> E*material-> Iy/le*(-6/le*d[2]-4*d[4]+6/le*d[8]-2*d[10]);
	stress[2]=material-> E*material-> Iy/le*(6/le*d[2]+2*d[4]-6/le*d[8]+4*d[10]);

	stress[3]=-material-> E*material-> Iy*(12/(le*le*le)*d[2]+6/(le*le)*d[4]
	          -12/(le*le*le)*d[8]+6/(le*le)*d[10]);

}
