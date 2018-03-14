/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Plate.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CPlate::CPlate()
{
	NEN = 4;	// Each element has 2 nodes
	nodes = new CNode*[NEN];
    
    ND = 12;
    LocationMatrix = new unsigned int[ND];

	ElementMaterial = nullptr;
}

//	Desconstructor
CPlate::~CPlate()
{
}

//	Read element data from stream Input
bool CPlate::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int N;

	Input >> N;	// element number

	if (N != Ele + 1)
	{
		cerr << "*** Error *** Elements must be inputted in order !" << endl 
			 << "    Expected element : " << Ele + 1 << endl
			 << "    Provided element : " << N << endl;

		return false;
	}

	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4;	// Four node numbers

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
    ElementMaterial = dynamic_cast<CPlateMaterial*>(MaterialSets) + MSet - 1;
	nodes[0] = &NodeList[N1 - 1];
	nodes[1] = &NodeList[N2 - 1];
	nodes[2] = &NodeList[N3 - 1];
	nodes[3] = &NodeList[N4 - 1];

	return true;
}

//	Write element data to stream
void CPlate::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele+1 << setw(11) << nodes[0]->NodeNumber
		   << setw(9) << nodes[1]->NodeNumber << setw(9) << nodes[2]->NodeNumber << setw(9) << nodes[3]->NodeNumber 
		   << setw(12) << ElementMaterial->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CPlate::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN; N++)
        for (unsigned int D = 2; D < 5; D++)
            LocationMatrix[i++] = nodes[N]->bcode[D];
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For Plate element, element stiffness is a 12x12 matrix, whose upper triangular part
//	has 78 elements
unsigned int CPlate::SizeOfStiffnessMatrix() { return 78; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CPlate::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());
	

// Elemenet material
	CPlateMaterial* material_ = (CPlateMaterial*)ElementMaterial;	// Pointer to material of the element

	double nu = material_->nu;
	double D0 = material_->E * material_->thick * material_->thick * material_->thick /12 / (1- nu * nu);
	double x[4],y[4],z[4],n[3],n2; // Coordinates of four nodes, they have to be the 4 nodes of a rectangle.
	for (unsigned int i = 0; i < 4; i++)
		{
			x[i]=nodes[i]->XYZ[0];
			y[i]=nodes[i]->XYZ[1];
			z[i]=nodes[i]->XYZ[2];
		}
	n[0] = (z[2] - z[0]) * (y[3] - y[1]) - (z[3] - z[1]) * (y[2] - y[0]);
	n[1] = (z[3] - z[1]) * (x[2] - x[0]) - (z[2] - z[0]) * (x[3] - x[1]);
	n[2] = (y[2] - y[0]) * (x[3] - x[1]) - (y[3] - y[1]) * (x[2] - x[0]);
	n2 = sqrt( n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	n[0] /= n2;
	n[1] /= n2;
	n[2] /= n2;
	n2 = (x[3] + x[1] - x[2] - x[0]) * n[0] + (y[3] + y[1] - y[2] - y[0]) * n[1] + (z[3] + z[1] - z[2] - z[0]) * n[2];
	double a = 0.25 * ( sqrt((x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0]) + (z[1] - z[0]) * (z[1] - z[0])) 
					  + sqrt((x[3] - x[2]) * (x[3] - x[2]) + (y[3] - y[2]) * (y[3] - y[2]) + (z[3] - z[2]) * (z[3] - z[2])) );
	double b = 0.25 * ( sqrt((x[1] - x[2]) * (x[1] - x[2]) + (y[1] - y[2]) * (y[1] - y[2]) + (z[1] - z[2]) * (z[1] - z[2])) 
					  + sqrt((x[3] - x[0]) * (x[3] - x[0]) + (y[3] - y[0]) * (y[3] - y[0]) + (z[3] - z[0]) * (z[3] - z[0])) );
	
	if ( abs(n2) > 1E-4 * a)
		cout << "Warning: element nodes are not in a plane" << endl;

	double a2 = a * a;
	double a3 = a2 * a;
	double a4 = a3 * a;
	double b2 = b * b;
	double b3 = b2 * b;
	double b4 = b3 * b;
	double ab = a * b;
	double a2b = a2 * b;
	double ab2 = a * b2;
	double a2b2 = a2 * b2;
	double a3b3 = a3 * b3;

	Matrix[0] = D0 * ( a4 + b4 + a2b2 * ( 7 - 2 * nu ) / 10 ) / a3b3;
	Matrix[1] = 4 * D0 * ( 5 * a2 - b2 * ( nu - 1 ) ) / 15 / ab;
	Matrix[2] = D0 * ( 1 + 10 * a2 / b2 + 4 * nu ) / 10 / a;
	Matrix[3] = 4 * D0 * ( 5 * b2 - a2 * ( nu -1) ) / 15 / ab;
	Matrix[4] = -D0 * nu;
	Matrix[5] = -D0 * ( b2 + a2 * ( 1 + 4 * nu ) / 10 ) / a2b;
	Matrix[6] = Matrix[0];
	Matrix[7] = D0 * ( 1 + 10 * b2 / a2 - nu ) / 10 / b;
	Matrix[8] = D0 * ( 5 * a2 - b2 * ( 1 + 4 * nu ) ) / 10 / ab2;
	Matrix[9] = D0 * ( 5 * a4 - 10 * b4 + a2b2 * ( 2 * nu - 7)) / 10 / a3b3;
	Matrix[10] = Matrix[1];
	Matrix[11] = Matrix[2];
	Matrix[12] = 0;
	Matrix[13] = 2 * D0 * ( 5 * a2 + 2 * b2 *( nu - 1 )) / 15 / ab;
	Matrix[14] = D0 * ( 5 * a2 - b2 * (1 + 4 * nu)) / 10 / ab2;
	Matrix[15] = Matrix[3];
	Matrix[16] = D0 * nu;
	Matrix[17] = D0 * ( 1 + 10 * b2 / a2 + 4 * nu) / 10 / b;
	Matrix[18] = D0 * (10 * b2 + a2 * ( nu - 1 )) / 15 / ab;
	Matrix[19] = 0;
	Matrix[20] = D0 * ( -10 * b2 + a2 * ( nu - 1 ) ) / 10 / a2b;
	Matrix[21] = Matrix[0];
	Matrix[22] = -D0 * ( -5 * b2 + a2 * ( 1 + 4 * nu)) / 10 / a2b;
	Matrix[23] = -D0 * ( 10 * a2 - b2 * ( nu - 1 ) ) / 10 / ab2;
	Matrix[24] = -D0 * ( 10 * a4 - 5 * b4 + a2b2 * (7 - 2 * nu ) ) / 10 /a3b3;
	Matrix[25] = D0 * ( 5 * b2 + a2 * ( nu - 1 ) ) / 10 / a2b;
	Matrix[26] = -D0 * (5 * a2 + b2 * ( nu - 1 ) ) / 10 /ab2;
	Matrix[27] = -D0 * (5 * a4 + 5 * b4 + a2b2 * ( 2 * nu - 7 ) ) / 10 / a3b3;
	Matrix[28] = Matrix[1];
	Matrix[29] = -Matrix[2];
	Matrix[30] = 0;
	Matrix[31] = D0 * ( 10 * a2 + b2 * ( nu - 1 ) ) / 15 / ab;
	Matrix[32] = D0 * ( 1 + 10 * a2 / b2 - nu ) / 10 / a;
	Matrix[33] = 0;
	Matrix[34] = D0 * ( 5 * a2 - b2 * ( nu - 1 ) ) / 15 / ab;
	Matrix[35] = D0 * ( 5 * a2 + b2 * ( nu - 1 ) ) / 10 / ab2;
	Matrix[36] = Matrix[3];
	Matrix[37] = -D0 * nu;
	Matrix[38] = D0 * ( 1 + 10 * b2 / a2 + 4 * nu ) / 10 / b;
	Matrix[39] = 2 * D0 * ( 5 * b2 + 2 * a2 * ( nu - 1 ) ) / 15 / ab;
	Matrix[40] = 0;
	Matrix[41] = -D0 * ( -5 * b2 + a2 * ( 1 + 4 * nu ) ) / 10 / a2b;
	Matrix[42] = D0 * ( 5 * b2 - a2 * ( nu - 1 ) ) / 15 / ab;
	Matrix[43] = 0;
	Matrix[44] = -D0 * ( 5 * b2 + a2 * ( nu - 1 ) ) / 10 / a2b;
	Matrix[45] = Matrix[0];
	Matrix[46] = D0 * ( -10 * b2 + a2 * ( nu - 1 ) ) / 10 / a2b;
	Matrix[47] = D0 * ( -5 * a2 + b2 * ( 1 + 4 * nu ) ) / 10 / ab2;
	Matrix[48] = D0 * ( 5 * a4 - 10 * b4 + a2b2 * ( 2 * nu - 7 ) ) / 10 / a3b3;
	Matrix[49] = -D0 * ( 5 * b2 + a2 * ( nu - 1 ) ) / 10 / a2b;
	Matrix[50] = -D0 * ( 5 * a2 + b2 * ( nu - 1 ) ) / 10 / ab2;
	Matrix[51] = -D0 * ( 5 * a4 + 5 * b4 + a2b2 * ( 2 * nu - 7 ) ) / 10 / a3b3;
	Matrix[52] = D0 * ( -5 * b2 + a2 * ( 1 + 4 * nu ) ) / 10 / a2b;
	Matrix[53] = -D0 * ( 10 * a2 - b2 * ( nu - 1 ) ) / 10 / ab2;
	Matrix[54] = -D0 * ( 10 * a4 - 5 * b4 + a2b2 * ( 7 - 2 * nu ) ) / 10 / a3b3;
	Matrix[55] = Matrix[1];
	Matrix[56] = -Matrix[2];
	Matrix[57] = 0;
	Matrix[58] = 2 * D0 * ( 5 * a2 + 2 * b2 * ( nu - 1 ) ) / 15 / ab;
	Matrix[59] = D0 * ( -5 * a2 + b2 * (1 + 4 * nu ) ) / 10 / ab2;
	Matrix[60] = 0;
	Matrix[61] = D0 * ( 5 * a2 - b2 * ( nu - 1 ) ) / 15 / ab;
	Matrix[62] = D0 * ( 5 * a2 + b2 * ( nu - 1 ) ) / 10 / ab2;
	Matrix[63] = 0;
	Matrix[64] = D0 * ( 10 * a2 + b2 * ( nu - 1 ) ) / 15 / ab;
	Matrix[65] = D0 * ( 1 + 10 * a2 / b2 - nu ) / 10 / a;
	Matrix[66] = Matrix[3];
	Matrix[67] = D0 * nu;
	Matrix[68] = -D0 * ( 10 * b2 + a2 * ( 1 + 4 * nu ) ) / 10 / a2b;
	Matrix[69] = D0 * ( 10 * b2 + a2 * ( nu - 1 ) ) / 15 / ab;
	Matrix[70] = 0;
	Matrix[71] = D0 * ( 1 + 10 * b2 / a2 - nu ) / 10 / b;
	Matrix[72] = D0 * ( 5 * b2 - a2 * ( nu - 1 ) ) / 15 / ab;
	Matrix[73] = 0;
	Matrix[74] = D0 * ( 5 * b2 + a2 * ( nu - 1 ) ) / 10 / a2b;
	Matrix[75] = 2 * D0 * ( 5 * b2 + 2 * a2 * ( nu - 1 ) ) / 15 / ab;
	Matrix[76] = 0;
	Matrix[77] = D0 * ( -5 * b2 + a2 * ( 1 + 4 * nu ) ) / 10 / a2b;
}

void CPlate::ElementGravity(double* bodyforce, double Gravity)
{
	clear(bodyforce, ND);
	

// Elemenet material
	CPlateMaterial* material_ = (CPlateMaterial*)ElementMaterial;	// Pointer to material of the element

	double x[4],y[4],z[4],n[3],n2; // Coordinates of four nodes, they have to be the 4 nodes of a rectangle.
	for (unsigned int i = 0; i < 4; i++)
		{
			x[i]=nodes[i]->XYZ[0];
			y[i]=nodes[i]->XYZ[1];
			z[i]=nodes[i]->XYZ[2];
		}
	n[0] = (z[2] - z[0]) * (y[3] - y[1]) - (z[3] - z[1]) * (y[2] - y[0]);
	n[1] = (z[3] - z[1]) * (x[2] - x[0]) - (z[2] - z[0]) * (x[3] - x[1]);
	n[2] = (y[2] - y[0]) * (x[3] - x[1]) - (y[3] - y[1]) * (x[2] - x[0]);
	n2 = sqrt( n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	n[0] /= n2;
	n[1] /= n2;
	n[2] /= n2;
	n2 = (x[3] + x[1] - x[2] - x[0]) * n[0] + (y[3] + y[1] - y[2] - y[0]) * n[1] + (z[3] + z[1] - z[2] - z[0]) * n[2];
	double a = 0.25 * ( sqrt((x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0]) + (z[1] - z[0]) * (z[1] - z[0])) 
					  + sqrt((x[3] - x[2]) * (x[3] - x[2]) + (y[3] - y[2]) * (y[3] - y[2]) + (z[3] - z[2]) * (z[3] - z[2])) );
	double b = 0.25 * ( sqrt((x[1] - x[2]) * (x[1] - x[2]) + (y[1] - y[2]) * (y[1] - y[2]) + (z[1] - z[2]) * (z[1] - z[2])) 
					  + sqrt((x[3] - x[0]) * (x[3] - x[0]) + (y[3] - y[0]) * (y[3] - y[0]) + (z[3] - z[0]) * (z[3] - z[0])) );
	
	double ab = a * b;
	double a2b = a * a * b;
	double ab2 = a * b * b;

	double q = material_->rho * Gravity * material_->thick;
	bodyforce[0] = ab * q;
	bodyforce[1] = ab2 * q / 3;
	bodyforce[2] = -a2b * q / 3;
	bodyforce[3] = ab * q;
	bodyforce[4] = ab2 * q / 3;
	bodyforce[5] = a2b * q / 3;
	bodyforce[6] = ab * q;
	bodyforce[7] = -ab2 * q / 3;
	bodyforce[8] = a2b * q / 3;
	bodyforce[9] = ab * q;
	bodyforce[10] = ab2 * q / 3;
	bodyforce[11] = -1/3 * a2b * q;
}

//	Calculate element stress 
void CPlate::ElementStress(double* stress, double* Displacement)
{
	CPlateMaterial* material_ = dynamic_cast<CPlateMaterial*>(ElementMaterial);	// Pointer to material of the element
	double E = material_->E;
	double nu = material_->nu;
	double t = material_->thick;
	double x[4],y[4],z[4]; // Coordinates of four nodes, they have to be the 4 nodes of a rectangle.
	for (unsigned int i = 0; i < 4; i++)
		{
			x[i]=nodes[i]->XYZ[0];
			y[i]=nodes[i]->XYZ[1];
			z[i]=nodes[i]->XYZ[2];
		}
	double a = 0.25 * ( sqrt((x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0]) + (z[1] - z[0]) * (z[1] - z[0])) 
					  + sqrt((x[3] - x[2]) * (x[3] - x[2]) + (y[3] - y[2]) * (y[3] - y[2]) + (z[3] - z[2]) * (z[3] - z[2])) );
	double b = 0.25 * ( sqrt((x[1] - x[2]) * (x[1] - x[2]) + (y[1] - y[2]) * (y[1] - y[2]) + (z[1] - z[2]) * (z[1] - z[2])) 
					  + sqrt((x[3] - x[0]) * (x[3] - x[0]) + (y[3] - y[0]) * (y[3] - y[0]) + (z[3] - z[0]) * (z[3] - z[0])) );
	
	double gauss[2] = {-1/sqrt(3),1/sqrt(3)};
	double xi,eta;
	double epsxx,epsyy,epsxy;

	unsigned int count = 0;
// Gauss integral, 2 points per dimension
	for (unsigned int i = 0; i < 2; i++)
	{
		for (unsigned int j = 0; j < 2; j++)
		{
			xi = gauss[i]; eta = gauss[j];

			double dis[12];
			for (unsigned int i = 0; i < 12; i++)
			{
				if (LocationMatrix[i])
					dis[i] = Displacement[LocationMatrix[i]-1];
				else
					dis[i] = 0;
			}
			epsxx = ( 3 * ( ( dis[3] - dis[0] ) * ( eta - 1 ) + ( dis[9] - dis[6] ) * ( 1 + eta ) ) * xi / ( 4 * a * a )
				    + ( ( eta - 1 ) * ( 3 * xi - 1 ) * dis[2] + ( eta - 1 ) * ( 3 * xi + 1 ) * dis[5]
					- ( eta + 1 ) * ( 3 * xi - 1 ) * dis[11] - ( eta + 1 ) * ( 3 * xi + 1 ) * dis[8] ) / 4 / a ) * t / 2;
			epsyy = ( 3 * eta * ( ( dis[9]- dis[0] ) * ( xi - 1 ) + ( dis[3] - dis[6] ) * ( xi + 1 ) ) / 4 / b / b 
					+ ( dis[1] * ( xi - 1 ) * ( 1 - 3 * eta ) + dis[4] * ( xi + 1 ) * ( 3 * eta -1 )
					+ dis[7] * ( xi + 1 ) * ( 3 * eta + 1 ) - dis[10] * ( xi - 1 ) * ( 3 * eta + 1 ) ) / 4 / b ) * t / 2;
			epsxy = ( ( ( 1 + 2 * eta - 3 * eta * eta ) * dis[1] - ( 1 + 2 * eta - 3 * eta * eta ) * dis[4] 
					+ ( -1 + 2 * eta + 3 * eta * eta ) * dis[7] + ( 1 - 2 * eta - 3 * eta * eta ) * dis[10] ) / 4 / a 
					+ (	dis[2] * ( -1 - 2 * xi + 3 * xi * xi ) + dis[5] * ( -1 + 2 * xi + 3 * xi * xi )
					+ dis[8] * ( 1 - 2 * xi - 3 * xi * xi ) + dis[11] * (1 + 2 * xi - 3 * xi * xi ) ) / 4 / b
					+( ( dis[0] + dis[6]) * (4 - 3 * eta * eta - 3 * xi * xi ) 
					+ ( dis[3] + dis[9] ) * ( -4 + 3 * eta * eta + 3 * xi * xi ) ) / 4 / a / b ) * t / 2;
			double k = E / (1 - nu * nu);
			stress[3*count] = k * (epsxx + nu * epsyy);
			stress[3*count+1] = k * (nu * epsxx + epsyy);
			stress[3*count+2] = k * (1 - nu) * epsxy / 2;
			count++;
		}
	}
}
