/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.0, October 14, 2017                                         */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/
/*                                                                           */
/*     Xiao Feiyu                                        */
/*     3T	2017                                                                      */
/*                                                    */
/*****************************************************************************/
#include "3T.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
C3T::C3T()
{
	NEN = 3;	// Each element has 4 nodes
	nodes = new CNode*[NEN];
    
    ND = 6;  //The Location Matrix is stored in one dimension
    LocationMatrix = new unsigned int[ND];

	ElementMaterial = NULL;

	location= NULL;
}

//	Desconstructor
C3T::~C3T()
{
	delete [] nodes;
    delete [] LocationMatrix;
}

//	Read element data from stream Input
bool C3T::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
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
	unsigned int N1, N2, N3;	// Left node number and right node number

	Input >> N1 >> N2 >> N3 >> MSet;
	ElementMaterial = &(dynamic_cast<C3TMaterial*>(MaterialSets))[MSet - 1];
	nodes[0] = &NodeList[N1 - 1];
	nodes[1] = &NodeList[N2 - 1];
	nodes[2] = &NodeList[N3 - 1];
	

	return true;
}

//	Write element data to stream OutputFile
void C3T::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele+1 << setw(9) << nodes[0]->NodeNumber
		 << setw(9) << nodes[1]->NodeNumber << setw(9) << nodes[2]->NodeNumber
		 << setw(12) << ElementMaterial->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void C3T::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN; N++)
        for (unsigned int D = 0; D < 2; D++)    // D=3 three dimensions
            LocationMatrix[i++] = nodes[N]->bcode[D];
}


//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 3 node 3T element, element stiffness is a 6x6 matrix, whose upper triangular part
//	has 21 elements
unsigned int C3T::SizeOfStiffnessMatrix() { return 21; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void C3T::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());
	//	Get the nodes --> C
	double C[3][2];
	for (unsigned int j=0; j<3; j++)
		for (unsigned int i = 0; i < 2; i++)
			C[j][i] = nodes[j]->XYZ[i] ;
	


//	Calculate element stiffness matrix

	C3TMaterial* material = dynamic_cast<C3TMaterial*>(ElementMaterial);	// Pointer to material of the element

	double E = material->E;
	double mu = material->mu;
    //  change here
    double scale;
    
    double Ke[6][6];
    // Get D
    //Plane stress
    double D_scale=E/(1-mu*mu);
    double D[3][3]={{D_scale,mu*D_scale,0},{mu*D_scale,
    	D_scale,0},{0,0,(1-mu)*D_scale/2}};
    //Plane strain
    /*double D_scale=E/((1-2*mu)*(1+mu));
    double D[3][3]={{D_scale*(1-mu),mu*D_scale,0},{mu*D_scale,
    	D_scale*(1-mu),0},{0,0,(1-2*mu)*D_scale/2}};*/
   
   	double B[3][6];
	
   	double Bt[6][3];
	double BtD[6][3];
   	double Kk[6][6];
	
   	double J[2][2];
   	double detJ;

   	double X1,X2,X3,Y1,Y2,Y3;
   	double B1x,B2x,B3x,B1y,B2y,B3y;
	
	X1=C[0][0];
	X2=C[1][0];
	X3=C[2][0];
	Y1=C[0][1];
	Y2=C[1][1];
	Y3=C[2][1];

	B1x=Y2-Y3;
	B2x=Y3-Y1;
	B3x=Y1-Y2;
	B1y=X3-X2;
	B2y=X1-X3;
	B3y=X2-X1;

	J[0][0]=X1-X3;
	J[0][1]=Y1-Y3;
	J[1][0]=X2-X3;
	J[1][1]=Y2-Y3;
  		
	detJ=J[0][0]*J[1][1]-J[0][1]*J[1][0];

	B[0][0]=B1x/detJ;
	B[0][1]=0;
	B[0][2]=B2x/detJ;
	B[0][3]=0;
	B[0][4]=B3x/detJ;
	B[0][5]=0;
	B[1][0]=0;
	B[1][1]=B1y/detJ;
	B[1][2]=0;
	B[1][3]=B2y/detJ;
	B[1][4]=0;
	B[1][5]=B3y/detJ;
	B[2][0]=B1y/detJ;
	B[2][1]=B1x/detJ;
	B[2][2]=B2y/detJ;
	B[2][3]=B2x/detJ;
	B[2][4]=B3y/detJ;
	B[2][5]=B3x/detJ;

	scale=0.5;

	for(unsigned m=0;m<6;m++)
				for (unsigned n=0;n<6;n++)
				{
					Ke[m][n]=0; 
				}

	//loop in gp
				//B[3][6]
				//transpose
				for(int ii=0;ii<3;ii++)
				{
					for(int jj=0;jj<6;jj++)
					{
						Bt[jj][ii]=B[ii][jj];
					}
				}
			//BtD[6][3]=Bt[6][3]*D[3][3]
		
			for(int ii=0;ii<6;++ii)
			{
			for(int jj=0;jj<3;++jj){
				BtD[ii][jj]=0;
					for(int kk=0;kk<3;++kk)
					{	
					BtD[ii][jj]+=Bt[ii][kk]*D[kk][jj];
					}
			}
			}
				//BtD[6][3]*B[3][6]
				
			for(int ii=0;ii<6;++ii)
			{
			for(int jj=0;jj<6;++jj){
				Kk[ii][jj]=0;
					for(int kk=0;kk<3;++kk)
					{
					Kk[ii][jj]+=BtD[ii][kk]*B[kk][jj];
					}
			}
			}
		
					
			for (int ii=0;ii<6;ii++)
				for(int jj=0;jj<6;jj++)
				{
					Ke[ii][jj]+=scale*Kk[ii][jj]*abs(detJ);
				}
		
	


	Matrix[0]=Ke[0][0];
	Matrix[1]=Ke[1][1];
	Matrix[2]=Ke[0][1];
	Matrix[3]=Ke[2][2];
	Matrix[4]=Ke[1][2];
	Matrix[5]=Ke[0][2];
	Matrix[6]=Ke[3][3];
	Matrix[7]=Ke[2][3];
	Matrix[8]=Ke[1][3];
	Matrix[9]=Ke[0][3];
	Matrix[10]=Ke[4][4];
	Matrix[11]=Ke[3][4];
	Matrix[12]=Ke[2][4];
	Matrix[13]=Ke[1][4];
	Matrix[14]=Ke[0][4];
	Matrix[15]=Ke[5][5];
	Matrix[16]=Ke[4][5];
	Matrix[17]=Ke[3][5];
	Matrix[18]=Ke[2][5];
	Matrix[19]=Ke[1][5];
	Matrix[20]=Ke[0][5];
}


void C3T::ElementGravity(double* bodyforce, double Gravity)
{
}

//	Calculate element stress 
void C3T::ElementStress(double* stress, double* Displacement)
{
	// ! The Guass points and weights
	int ngp=3 ; // Guass points of 4Q elements in each dimension
	// Guass point of 3
	double gp1[3]={0.16666666666,0.66666666666,0.16666666666};
	double gp2[3]={0.16666666666,0.16666666666,0.66666666666};
	// Gauss weight 
	//double w[3]={0.16666666666,0.16666666666,0.16666666666};
	

	//	Get the nodes --> C
	double C[3][2];
	for (unsigned int j=0; j<3; j++)
		for (unsigned int i = 0; i < 2; i++)
			C[j][i] = nodes[j]->XYZ[i] ;
	


//	Calculate element stiffness matrix

	C3TMaterial* material = dynamic_cast<C3TMaterial*>(ElementMaterial);	// Pointer to material of the element

	double E = material->E;
	double mu = material->mu;
    //  change here
    double eta;
    double psi;
  
  
    // Get D
    //Plane stress
    double D_scale=E/(1-mu*mu);
    double D[3][3]={{D_scale,mu*D_scale,0},{mu*D_scale,
    	D_scale,0},{0,0,(1-mu)*D_scale/2}};
    //Plane strain
    /*double D_scale=E/((1-2*mu)*(1+mu));
    double D[3][3]={{D_scale*(1-mu),mu*D_scale,0},{mu*D_scale,
    	D_scale*(1-mu),0},{0,0,(1-2*mu)*D_scale/2}};*/
   
  
    double B[3][6];
    double Na[3];
	double XY[2]={0,0};
    double s[3];//strain
	double S[3];//stress
	
	double J[2][2];
	double detJ;
	double N1,N2,N3;
	
	
	int it_times=0;

	double X1,X2,X3,Y1,Y2,Y3;
   	double B1x,B2x,B3x,B1y,B2y,B3y;
	
	X1=C[0][0];
	X2=C[1][0];
	X3=C[2][0];
	Y1=C[0][1];
	Y2=C[1][1];
	Y3=C[2][1];

	B1x=Y2-Y3;
	B2x=Y3-Y1;
	B3x=Y1-Y2;
	B1y=X3-X2;
	B2y=X1-X3;
	B3y=X2-X1;

	J[0][0]=X1-X3;
	J[0][1]=Y1-Y3;
	J[1][0]=X2-X3;
	J[1][1]=Y2-Y3;
  		
	detJ=J[0][0]*J[1][1]-J[0][1]*J[1][0];

	B[0][0]=B1x/detJ;
	B[0][1]=0;
	B[0][2]=B2x/detJ;
	B[0][3]=0;
	B[0][4]=B3x/detJ;
	B[0][5]=0/detJ;
	B[1][0]=0;
	B[1][1]=B1y/detJ;
	B[1][2]=0;
	B[1][3]=B2y/detJ;
	B[1][4]=0;
	B[1][5]=B3y/detJ;
	B[2][0]=B1y/detJ;
	B[2][1]=B1x/detJ;
	B[2][2]=B2y/detJ;
	B[2][3]=B2x/detJ;
	B[2][4]=B3y/detJ;
	B[2][5]=B3x/detJ;

	double x1,x2,x3,y1,y2,y3;
	x1=1.0;
	x2=0.0;
	x3=0.0;
	y1=0.0;
	y2=1.0;
	y3=0.0;

	// compute strains and stress at the gauss points  
    for(int i=0;i<ngp;i++)
    	{
			it_times+=1;
    		eta=gp1[i];
    		psi=gp2[i];


			N1 =  (x2*y3-x3*y2)+(y2-y3)*psi+(x3-x2)*eta;
			N2 =  (x3*y1-x1*y3)+(y3-y1)*psi+(x1-x3)*eta; 
			N3 =  (x1*y2-x2*y1)+(y1-y2)*psi+(x2-x1)*eta;
			

		
    		Na[0]=N1;
			Na[1]=N2;
			Na[2]=N3;
			
			
			//XY=Na[1][3]*C[3][2]
			XY[0]=0;
			XY[1]=0;
			
			for(int jj=0;jj<2;++jj){
				XY[jj]=0;
					for(int kk=0;kk<3;++kk)
					{
					XY[jj]+=Na[kk]*C[kk][jj];
					}
					cout<<"XY"<<XY[jj]<<endl;
			}
			//*(X_Guass+2*it_times-2)=XY[1];
			//*(X_Guass+2*it_times-1)=XY[2];
			
    					
				double de[6];

				for (int ii=0;ii<6;ii++)
				{
					de[ii]=0;
					if(LocationMatrix[ii])
					de[ii]=*(Displacement+(LocationMatrix[ii]-1));
				}
				//strain
				//s=B[3][6]*Dispalcement[6]
				for(int ii=0;ii<3;++ii)
					{
						s[ii]=0;
							for(int kk=0;kk<6;++kk)
							{	
							s[ii]+=B[ii][kk]*de[kk];			
							}
							//cout<<"s"<<s[ii]<<endl;
					}
				//stress
				//S[3]=D[3][3]*s[3]
				for(int ii=0;ii<3;++ii)
					{
						S[ii]=0;
							for(int kk=0;kk<3;++kk)
							{	
							S[ii]+=D[ii][kk]*s[kk];		
							}
					}
					
				for (int ii=3*(it_times-1);ii<3*it_times;ii++)
					*(stress+ii) = 0.0;
					

    		for (int ii = 0; ii < 3; ii++)
				{
					//if (LocationMatrix[ii])
						*(stress+3*(it_times-1)+ii) = S[ii] ;
				}

			for (int ii = 0; ii < 2; ii++)
			           *(location+2*(it_times-1)+ii) = XY[ii];
	
    	}	

		
}




        
       
        
