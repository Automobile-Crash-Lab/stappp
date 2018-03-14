/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "8H.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
C8H::C8H()
{
	NEN = 8;
	nodes = new CNode*[NEN];
    
    ND = 24;
    LocationMatrix = new unsigned int[ND];

	ElementMaterial = nullptr;
}

//	Desconstructor
C8H::~C8H()
{
	delete [] nodes;
    delete [] LocationMatrix;
}
//	Read element data from stream Input
bool C8H::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
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

	unsigned int MSet;
	unsigned int N1, N2, N3 ,N4 ,N5 ,N6 ,N7 ,N8;

	Input >> N1 >> N2 >> N3 >> N4 >> N5 >> N6 >> N7 >> N8 >> MSet;
    ElementMaterial = dynamic_cast<C8HMaterial*>(MaterialSets) + MSet - 1;
	nodes[0] = &NodeList[N1 - 1];
	nodes[1] = &NodeList[N2 - 1];
	nodes[2] = &NodeList[N3 - 1];
	nodes[3] = &NodeList[N4 - 1];
	nodes[4] = &NodeList[N5 - 1];
	nodes[5] = &NodeList[N6 - 1];
	nodes[6] = &NodeList[N7 - 1];
	nodes[7] = &NodeList[N8 - 1];

	return true;
}

//	Write element data to stream
void C8H::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele+1 
		 << setw(11) << nodes[0]->NodeNumber 
		 << setw(9) << nodes[1]->NodeNumber 
		 << setw(9) << nodes[2]->NodeNumber 
		 << setw(9) << nodes[3]->NodeNumber
		 << setw(9) << nodes[4]->NodeNumber
		 << setw(9) << nodes[5]->NodeNumber
		 << setw(9) << nodes[6]->NodeNumber
		 << setw(9) << nodes[7]->NodeNumber
		 << setw(12) << ElementMaterial->nset 
		 << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void C8H::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN; N++)
        for (unsigned int D = 0; D < 3; D++)
            LocationMatrix[i++] = nodes[N]->bcode[D];
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For Plate element, element stiffness is a 12x12 matrix, whose upper triangular part
//	has 78 elements
unsigned int C8H::SizeOfStiffnessMatrix() { return 300; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void C8H::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	C8HMaterial* material = dynamic_cast<C8HMaterial*>(ElementMaterial);
	double E8H = material->E;
	double nu8H = material->posi_ratio;
	double G0 = E8H/(1+nu8H);
	double lambda = nu8H*E8H/((1+nu8H)*(1-2*nu8H));

	double GPoint[2];
	GPoint[0] = -1/sqrt(3);
	GPoint[1] = 1/sqrt(3);
	for (unsigned int m = 0; m < 2; m++)
	{
		for (unsigned int n = 0; n < 2; n++)
		{
			for (unsigned int o = 0; o < 2; o++)
			{
				double xi = GPoint[m];
				double eta = GPoint[n];
				double zeta = GPoint[o];
				double GN[12];
				GN[0] = 0.125*(1 - eta)*(1 - zeta);
				GN[1] = 0.125*(1 + eta)*(1 - zeta);
				GN[2] = 0.125*(1 - eta)*(1 + zeta);
				GN[3] = 0.125*(1 + eta)*(1 + zeta);
				GN[4] = 0.125*(1 - xi)*(1 - zeta);
				GN[5] = 0.125*(1 + xi)*(1 - zeta);
				GN[6] = 0.125*(1 - xi)*(1 + zeta);
				GN[7] = 0.125*(1 + xi)*(1 + zeta);
				GN[8] = 0.125*(1 - xi)*(1 - eta);
				GN[9] = 0.125*(1 + xi)*(1 - eta);
				GN[10] = 0.125*(1 + xi)*(1 + eta);
				GN[11] = 0.125*(1 - xi)*(1 + eta);

				double J[9];
				J[0] = GN[0]*(nodes[1]->XYZ[0] - nodes[0]->XYZ[0]) + GN[1]*(nodes[2]->XYZ[0] - nodes[3]->XYZ[0]) + GN[2]*(nodes[5]->XYZ[0] - nodes[4]->XYZ[0]) + GN[3]*(nodes[6]->XYZ[0] - nodes[7]->XYZ[0]);
				J[1] = GN[0]*(nodes[1]->XYZ[1] - nodes[0]->XYZ[1]) + GN[1]*(nodes[2]->XYZ[1] - nodes[3]->XYZ[1]) + GN[2]*(nodes[5]->XYZ[1] - nodes[4]->XYZ[1]) + GN[3]*(nodes[6]->XYZ[1] - nodes[7]->XYZ[1]);
				J[2] = GN[0]*(nodes[1]->XYZ[2] - nodes[0]->XYZ[2]) + GN[1]*(nodes[2]->XYZ[2] - nodes[3]->XYZ[2]) + GN[2]*(nodes[5]->XYZ[2] - nodes[4]->XYZ[2]) + GN[3]*(nodes[6]->XYZ[2] - nodes[7]->XYZ[2]);
				J[3] = GN[4]*(nodes[3]->XYZ[0] - nodes[0]->XYZ[0]) + GN[5]*(nodes[2]->XYZ[0] - nodes[1]->XYZ[0]) + GN[6]*(nodes[7]->XYZ[0] - nodes[4]->XYZ[0]) + GN[7]*(nodes[6]->XYZ[0] - nodes[5]->XYZ[0]);
				J[4] = GN[4]*(nodes[3]->XYZ[1] - nodes[0]->XYZ[1]) + GN[5]*(nodes[2]->XYZ[1] - nodes[1]->XYZ[1]) + GN[6]*(nodes[7]->XYZ[1] - nodes[4]->XYZ[1]) + GN[7]*(nodes[6]->XYZ[1] - nodes[5]->XYZ[1]);
				J[5] = GN[4]*(nodes[3]->XYZ[2] - nodes[0]->XYZ[2]) + GN[5]*(nodes[2]->XYZ[2] - nodes[1]->XYZ[2]) + GN[6]*(nodes[7]->XYZ[2] - nodes[4]->XYZ[2]) + GN[7]*(nodes[6]->XYZ[2] - nodes[5]->XYZ[2]);
				J[6] = GN[8]*(nodes[4]->XYZ[0] - nodes[0]->XYZ[0]) + GN[9]*(nodes[5]->XYZ[0] - nodes[1]->XYZ[0]) + GN[10]*(nodes[6]->XYZ[0] - nodes[2]->XYZ[0]) + GN[11]*(nodes[7]->XYZ[0] - nodes[3]->XYZ[0]);
				J[7] = GN[8]*(nodes[4]->XYZ[1] - nodes[0]->XYZ[1]) + GN[9]*(nodes[5]->XYZ[1] - nodes[1]->XYZ[1]) + GN[10]*(nodes[6]->XYZ[1] - nodes[2]->XYZ[1]) + GN[11]*(nodes[7]->XYZ[1] - nodes[3]->XYZ[1]);
				J[8] = GN[8]*(nodes[4]->XYZ[2] - nodes[0]->XYZ[2]) + GN[9]*(nodes[5]->XYZ[2] - nodes[1]->XYZ[2]) + GN[10]*(nodes[6]->XYZ[2] - nodes[2]->XYZ[2]) + GN[11]*(nodes[7]->XYZ[2] - nodes[3]->XYZ[2]);
				double DETJ1 = J[0]*J[4]*J[8] - J[0]*J[5]*J[7] - J[1]*J[3]*J[8] + J[1]*J[5]*J[6] + J[2]*J[3]*J[7] - J[2]*J[4]*J[6];
				double DETJ = abs(DETJ1);

				double INVJ[9];
				INVJ[0] = (J[4]*J[8]-J[5]*J[7])/DETJ1;
			    INVJ[1] = -(J[1]*J[8]-J[2]*J[7])/DETJ1;
			    INVJ[2] = (J[1]*J[5]-J[2]*J[4])/DETJ1;
			    INVJ[3] = -(J[3]*J[8]-J[5]*J[6])/DETJ1;
			    INVJ[4] = (J[0]*J[8]-J[2]*J[6])/DETJ1;
			    INVJ[5] = -(J[0]*J[5]-J[2]*J[3])/DETJ1;
			    INVJ[6] = (J[3]*J[7]-J[4]*J[6])/DETJ1;
			    INVJ[7] = -(J[0]*J[7]-J[1]*J[6])/DETJ1;
			    INVJ[8] = (J[0]*J[4]-J[1]*J[3])/DETJ1;

				double B1[24];
				B1[0] = -GN[0]*INVJ[0]-GN[4]*INVJ[1]-GN[8]*INVJ[2];
			    B1[1] = GN[0]*INVJ[0]-GN[5]*INVJ[1]-GN[9]*INVJ[2];
			    B1[2] = GN[1]*INVJ[0]+GN[5]*INVJ[1]-GN[10]*INVJ[2];
			    B1[3] = -GN[1]*INVJ[0]+GN[4]*INVJ[1]-GN[11]*INVJ[2];
				B1[4] = -GN[2]*INVJ[0]-GN[6]*INVJ[1]+GN[8]*INVJ[2];
				B1[5] = GN[2]*INVJ[0]-GN[7]*INVJ[1]+GN[9]*INVJ[2];
				B1[6] = GN[3]*INVJ[0]+GN[7]*INVJ[1]+GN[10]*INVJ[2];
				B1[7] = -GN[3]*INVJ[0]+GN[6]*INVJ[1]+GN[11]*INVJ[2];
			    B1[8] = -GN[0]*INVJ[3]-GN[4]*INVJ[4]-GN[8]*INVJ[5];
			    B1[9] = GN[0]*INVJ[3]-GN[5]*INVJ[4]-GN[9]*INVJ[5];
			    B1[10] = GN[1]*INVJ[3]+GN[5]*INVJ[4]-GN[10]*INVJ[5];
			    B1[11] = -GN[1]*INVJ[3]+GN[4]*INVJ[4]-GN[11]*INVJ[5];
				B1[12] = -GN[2]*INVJ[3]-GN[6]*INVJ[4]+GN[8]*INVJ[5];
				B1[13] = GN[2]*INVJ[3]-GN[7]*INVJ[4]+GN[9]*INVJ[5];
				B1[14] = GN[3]*INVJ[3]+GN[7]*INVJ[4]+GN[10]*INVJ[5];
				B1[15] = -GN[3]*INVJ[3]+GN[6]*INVJ[4]+GN[11]*INVJ[5];
				B1[16] = -GN[0]*INVJ[6]-GN[4]*INVJ[7]-GN[8]*INVJ[8];
			    B1[17] = GN[0]*INVJ[6]-GN[5]*INVJ[7]-GN[9]*INVJ[8];
			    B1[18] = GN[1]*INVJ[6]+GN[5]*INVJ[7]-GN[10]*INVJ[8];
			    B1[19] = -GN[1]*INVJ[6]+GN[4]*INVJ[7]-GN[11]*INVJ[8];
				B1[20] = -GN[2]*INVJ[6]-GN[6]*INVJ[7]+GN[8]*INVJ[8];
				B1[21] = GN[2]*INVJ[6]-GN[7]*INVJ[7]+GN[9]*INVJ[8];
				B1[22] = GN[3]*INVJ[6]+GN[7]*INVJ[7]+GN[10]*INVJ[8];
				B1[23] = -GN[3]*INVJ[6]+GN[6]*INVJ[7]+GN[11]*INVJ[8];


			       Matrix[0] += (B1[0]*B1[0]*(G0+lambda)+B1[8]*G0*B1[8]+B1[16]*G0*B1[16])*DETJ;
				   Matrix[1] += (B1[8]*B1[8]*(G0+lambda)+B1[0]*G0*B1[0]+B1[16]*G0*B1[16])*DETJ;
				   Matrix[2] += (B1[0]*G0*B1[8]+B1[8]*lambda*B1[0])*DETJ;
				   Matrix[3] += (B1[16]*B1[16]*(G0+lambda)+B1[0]*G0*B1[0]+B1[8]*G0*B1[8])*DETJ;
				   Matrix[4] += (B1[8]*G0*B1[16]+B1[16]*lambda*B1[8])*DETJ;
				   Matrix[5] += (B1[0]*G0*B1[16]+B1[16]*lambda*B1[0])*DETJ;
				   Matrix[6] += (B1[1]*B1[1]*(G0+lambda)+B1[9]*G0*B1[9]+B1[17]*G0*B1[17])*DETJ;
				   Matrix[7] += (B1[17]*G0*B1[0]+B1[1]*lambda*B1[16])*DETJ;
				   Matrix[8] += (B1[9]*G0*B1[0]+B1[1]*lambda*B1[8])*DETJ;
				   Matrix[9] += (B1[1]*B1[0]*(G0+lambda)+B1[9]*G0*B1[8]+B1[17]*G0*B1[16])*DETJ;
				   Matrix[10] += (B1[9]*B1[9]*(G0+lambda)+B1[1]*G0*B1[1]+B1[17]*G0*B1[17])*DETJ;
				   Matrix[11] += (B1[1]*G0*B1[9]+B1[9]*lambda*B1[1])*DETJ;
				   Matrix[12] += (B1[17]*G0*B1[8]+B1[9]*lambda*B1[16])*DETJ;
				   Matrix[13] += (B1[9]*B1[8]*(G0+lambda)+B1[1]*G0*B1[0]+B1[17]*G0*B1[16])*DETJ;
				   Matrix[14] += (B1[1]*G0*B1[8]+B1[9]*lambda*B1[0])*DETJ;
				   Matrix[15] += (B1[17]*B1[17]*(G0+lambda)+B1[1]*G0*B1[1]+B1[9]*G0*B1[9])*DETJ;
				   Matrix[16] += (B1[9]*G0*B1[17]+B1[17]*lambda*B1[9])*DETJ;
				   Matrix[17] += (B1[1]*G0*B1[17]+B1[17]*lambda*B1[1])*DETJ;
				   Matrix[18] += (B1[17]*B1[16]*(G0+lambda)+B1[1]*G0*B1[0]+B1[9]*G0*B1[8])*DETJ;
				   Matrix[19] += (B1[9]*G0*B1[16]+B1[17]*lambda*B1[8])*DETJ;
				   Matrix[20] += (B1[1]*G0*B1[16]+B1[17]*lambda*B1[0])*DETJ;
				   Matrix[21] += (B1[2]*B1[2]*(G0+lambda)+B1[10]*G0*B1[10]+B1[18]*G0*B1[18])*DETJ;
				   Matrix[22] += (B1[18]*G0*B1[1]+B1[2]*lambda*B1[17])*DETJ;
				   Matrix[23] += (B1[10]*G0*B1[1]+B1[2]*lambda*B1[9])*DETJ;
				   Matrix[24] += (B1[2]*B1[1]*(G0+lambda)+B1[10]*G0*B1[9]+B1[18]*G0*B1[17])*DETJ;
				   Matrix[25] += (B1[18]*G0*B1[0]+B1[2]*lambda*B1[16])*DETJ;
				   Matrix[26] += (B1[10]*G0*B1[0]+B1[2]*lambda*B1[8])*DETJ;
				   Matrix[27] += (B1[2]*B1[0]*(G0+lambda)+B1[10]*G0*B1[8]+B1[18]*G0*B1[16])*DETJ;
				   Matrix[28] += (B1[10]*B1[10]*(G0+lambda)+B1[2]*G0*B1[2]+B1[18]*G0*B1[18])*DETJ;
				   Matrix[29] += (B1[2]*G0*B1[10]+B1[10]*lambda*B1[2])*DETJ;
				   Matrix[30] += (B1[18]*G0*B1[9]+B1[10]*lambda*B1[17])*DETJ;
				   Matrix[31] += (B1[10]*B1[9]*(G0+lambda)+B1[2]*G0*B1[1]+B1[18]*G0*B1[17])*DETJ;
				   Matrix[32] += (B1[2]*G0*B1[9]+B1[10]*lambda*B1[1])*DETJ;
				   Matrix[33] += (B1[18]*G0*B1[8]+B1[10]*lambda*B1[16])*DETJ;
				   Matrix[34] += (B1[10]*B1[8]*(G0+lambda)+B1[2]*G0*B1[0]+B1[18]*G0*B1[16])*DETJ;
				   Matrix[35] += (B1[2]*G0*B1[8]+B1[10]*lambda*B1[0])*DETJ;
				   Matrix[36] += (B1[18]*B1[18]*(G0+lambda)+B1[2]*G0*B1[2]+B1[10]*G0*B1[10])*DETJ;
				   Matrix[37] += (B1[10]*G0*B1[18]+B1[18]*lambda*B1[10])*DETJ;
				   Matrix[38] += (B1[2]*G0*B1[18]+B1[18]*lambda*B1[2])*DETJ;
				   Matrix[39] += (B1[18]*B1[17]*(G0+lambda)+B1[2]*G0*B1[1]+B1[10]*G0*B1[9])*DETJ;
				   Matrix[40] += (B1[10]*G0*B1[17]+B1[18]*lambda*B1[9])*DETJ;
				   Matrix[41] += (B1[2]*G0*B1[17]+B1[18]*lambda*B1[1])*DETJ;
				   Matrix[42] += (B1[18]*B1[16]*(G0+lambda)+B1[2]*G0*B1[0]+B1[10]*G0*B1[8])*DETJ;
				   Matrix[43] += (B1[10]*G0*B1[16]+B1[18]*lambda*B1[8])*DETJ;
				   Matrix[44] += (B1[2]*G0*B1[16]+B1[18]*lambda*B1[0])*DETJ;
				   Matrix[45] += (B1[3]*B1[3]*(G0+lambda)+B1[11]*G0*B1[11]+B1[19]*G0*B1[19])*DETJ;
				   Matrix[46] += (B1[19]*G0*B1[2]+B1[3]*lambda*B1[18])*DETJ;
				   Matrix[47] += (B1[11]*G0*B1[2]+B1[3]*lambda*B1[10])*DETJ;
				   Matrix[48] += (B1[3]*B1[2]*(G0+lambda)+B1[11]*G0*B1[10]+B1[19]*G0*B1[18])*DETJ;
				   Matrix[49] += (B1[19]*G0*B1[1]+B1[3]*lambda*B1[17])*DETJ;
				   Matrix[50] += (B1[11]*G0*B1[1]+B1[3]*lambda*B1[9])*DETJ;
				   Matrix[51] += (B1[3]*B1[1]*(G0+lambda)+B1[11]*G0*B1[9]+B1[19]*G0*B1[17])*DETJ;
				   Matrix[52] += (B1[19]*G0*B1[0]+B1[3]*lambda*B1[16])*DETJ;
				   Matrix[53] += (B1[11]*G0*B1[0]+B1[3]*lambda*B1[8])*DETJ;
				   Matrix[54] += (B1[3]*B1[0]*(G0+lambda)+B1[11]*G0*B1[8]+B1[19]*G0*B1[16])*DETJ;
				   Matrix[55] += (B1[11]*B1[11]*(G0+lambda)+B1[3]*G0*B1[3]+B1[19]*G0*B1[19])*DETJ;
				   Matrix[56] += (B1[3]*G0*B1[11]+B1[11]*lambda*B1[3])*DETJ;
				   Matrix[57] += (B1[19]*G0*B1[10]+B1[11]*lambda*B1[18])*DETJ;
				   Matrix[58] += (B1[11]*B1[10]*(G0+lambda)+B1[3]*G0*B1[2]+B1[19]*G0*B1[18])*DETJ;
				   Matrix[59] += (B1[3]*G0*B1[10]+B1[11]*lambda*B1[2])*DETJ;
				   Matrix[60] += (B1[19]*G0*B1[9]+B1[11]*lambda*B1[17])*DETJ;
				   Matrix[61] += (B1[11]*B1[9]*(G0+lambda)+B1[3]*G0*B1[1]+B1[19]*G0*B1[17])*DETJ;
				   Matrix[62] += (B1[3]*G0*B1[9]+B1[11]*lambda*B1[1])*DETJ;
				   Matrix[63] += (B1[19]*G0*B1[8]+B1[11]*lambda*B1[16])*DETJ;
				   Matrix[64] += (B1[11]*B1[8]*(G0+lambda)+B1[3]*G0*B1[0]+B1[19]*G0*B1[16])*DETJ;
				   Matrix[65] += (B1[3]*G0*B1[8]+B1[11]*lambda*B1[0])*DETJ;
				   Matrix[66] += (B1[19]*B1[19]*(G0+lambda)+B1[3]*G0*B1[3]+B1[11]*G0*B1[11])*DETJ;
				   Matrix[67] += (B1[11]*G0*B1[19]+B1[19]*lambda*B1[11])*DETJ;
				   Matrix[68] += (B1[3]*G0*B1[19]+B1[19]*lambda*B1[3])*DETJ;
				   Matrix[69] += (B1[19]*B1[18]*(G0+lambda)+B1[3]*G0*B1[2]+B1[11]*G0*B1[10])*DETJ;
				   Matrix[70] += (B1[11]*G0*B1[18]+B1[19]*lambda*B1[10])*DETJ;
				   Matrix[71] += (B1[3]*G0*B1[18]+B1[19]*lambda*B1[2])*DETJ;
				   Matrix[72] += (B1[19]*B1[17]*(G0+lambda)+B1[3]*G0*B1[1]+B1[11]*G0*B1[9])*DETJ;
				   Matrix[73] += (B1[11]*G0*B1[17]+B1[19]*lambda*B1[9])*DETJ;
				   Matrix[74] += (B1[3]*G0*B1[17]+B1[19]*lambda*B1[1])*DETJ;
				   Matrix[75] += (B1[19]*B1[16]*(G0+lambda)+B1[3]*G0*B1[0]+B1[11]*G0*B1[8])*DETJ;
				   Matrix[76] += (B1[11]*G0*B1[16]+B1[19]*lambda*B1[8])*DETJ;
				   Matrix[77] += (B1[3]*G0*B1[16]+B1[19]*lambda*B1[0])*DETJ;
				   Matrix[78] += (B1[4]*B1[4]*(G0+lambda)+B1[12]*G0*B1[12]+B1[20]*G0*B1[20])*DETJ;
				   Matrix[79] += (B1[20]*G0*B1[3]+B1[4]*lambda*B1[19])*DETJ;
				   Matrix[80] += (B1[12]*G0*B1[3]+B1[4]*lambda*B1[11])*DETJ;
				   Matrix[81] += (B1[4]*B1[3]*(G0+lambda)+B1[12]*G0*B1[11]+B1[20]*G0*B1[19])*DETJ;
				   Matrix[82] += (B1[20]*G0*B1[2]+B1[4]*lambda*B1[18])*DETJ;
				   Matrix[83] += (B1[12]*G0*B1[2]+B1[4]*lambda*B1[10])*DETJ;
				   Matrix[84] += (B1[4]*B1[2]*(G0+lambda)+B1[12]*G0*B1[10]+B1[20]*G0*B1[18])*DETJ;
				   Matrix[85] += (B1[20]*G0*B1[1]+B1[4]*lambda*B1[17])*DETJ;
				   Matrix[86] += (B1[12]*G0*B1[1]+B1[4]*lambda*B1[9])*DETJ;
				   Matrix[87] += (B1[4]*B1[1]*(G0+lambda)+B1[12]*G0*B1[9]+B1[20]*G0*B1[17])*DETJ;
				   Matrix[88] += (B1[20]*G0*B1[0]+B1[4]*lambda*B1[16])*DETJ;
				   Matrix[89] += (B1[12]*G0*B1[0]+B1[4]*lambda*B1[8])*DETJ;
				   Matrix[90] += (B1[4]*B1[0]*(G0+lambda)+B1[12]*G0*B1[8]+B1[20]*G0*B1[16])*DETJ;
				   Matrix[91] += (B1[12]*B1[12]*(G0+lambda)+B1[4]*G0*B1[4]+B1[20]*G0*B1[20])*DETJ;
				   Matrix[92] += (B1[4]*G0*B1[12]+B1[12]*lambda*B1[4])*DETJ;
				   Matrix[93] += (B1[20]*G0*B1[11]+B1[12]*lambda*B1[19])*DETJ;
				   Matrix[94] += (B1[12]*B1[11]*(G0+lambda)+B1[4]*G0*B1[3]+B1[20]*G0*B1[19])*DETJ;
				   Matrix[95] += (B1[4]*G0*B1[11]+B1[12]*lambda*B1[3])*DETJ;
				   Matrix[96] += (B1[20]*G0*B1[10]+B1[12]*lambda*B1[18])*DETJ;
				   Matrix[97] += (B1[12]*B1[10]*(G0+lambda)+B1[4]*G0*B1[2]+B1[20]*G0*B1[18])*DETJ;
				   Matrix[98] += (B1[4]*G0*B1[10]+B1[12]*lambda*B1[2])*DETJ;
				   Matrix[99] += (B1[20]*G0*B1[9]+B1[12]*lambda*B1[17])*DETJ;
				   Matrix[100] += (B1[12]*B1[9]*(G0+lambda)+B1[4]*G0*B1[1]+B1[20]*G0*B1[17])*DETJ;
				   Matrix[101] += (B1[4]*G0*B1[9]+B1[12]*lambda*B1[1])*DETJ;
				   Matrix[102] += (B1[20]*G0*B1[8]+B1[12]*lambda*B1[16])*DETJ;
				   Matrix[103] += (B1[12]*B1[8]*(G0+lambda)+B1[4]*G0*B1[0]+B1[20]*G0*B1[16])*DETJ;
				   Matrix[104] += (B1[4]*G0*B1[8]+B1[12]*lambda*B1[0])*DETJ;
				   Matrix[105] += (B1[20]*B1[20]*(G0+lambda)+B1[4]*G0*B1[4]+B1[12]*G0*B1[12])*DETJ;
				   Matrix[106] += (B1[12]*G0*B1[20]+B1[20]*lambda*B1[12])*DETJ;
				   Matrix[107] += (B1[4]*G0*B1[20]+B1[20]*lambda*B1[4])*DETJ;
				   Matrix[108] += (B1[20]*B1[19]*(G0+lambda)+B1[4]*G0*B1[3]+B1[12]*G0*B1[11])*DETJ;
				   Matrix[109] += (B1[12]*G0*B1[19]+B1[20]*lambda*B1[11])*DETJ;
				   Matrix[110] += (B1[4]*G0*B1[19]+B1[20]*lambda*B1[3])*DETJ;
				   Matrix[111] += (B1[20]*B1[18]*(G0+lambda)+B1[4]*G0*B1[2]+B1[12]*G0*B1[10])*DETJ;
				   Matrix[112] += (B1[12]*G0*B1[18]+B1[20]*lambda*B1[10])*DETJ;
				   Matrix[113] += (B1[4]*G0*B1[18]+B1[20]*lambda*B1[2])*DETJ;
				   Matrix[114] += (B1[20]*B1[17]*(G0+lambda)+B1[4]*G0*B1[1]+B1[12]*G0*B1[9])*DETJ;
				   Matrix[115] += (B1[12]*G0*B1[17]+B1[20]*lambda*B1[9])*DETJ;
				   Matrix[116] += (B1[4]*G0*B1[17]+B1[20]*lambda*B1[1])*DETJ;
				   Matrix[117] += (B1[20]*B1[16]*(G0+lambda)+B1[4]*G0*B1[0]+B1[12]*G0*B1[8])*DETJ;
				   Matrix[118] += (B1[12]*G0*B1[16]+B1[20]*lambda*B1[8])*DETJ;
				   Matrix[119] += (B1[4]*G0*B1[16]+B1[20]*lambda*B1[0])*DETJ;
				   Matrix[120] += (B1[5]*B1[5]*(G0+lambda)+B1[13]*G0*B1[13]+B1[21]*G0*B1[21])*DETJ;
				   Matrix[121] += (B1[21]*G0*B1[4]+B1[5]*lambda*B1[20])*DETJ;
				   Matrix[122] += (B1[13]*G0*B1[4]+B1[5]*lambda*B1[12])*DETJ;
				   Matrix[123] += (B1[5]*B1[4]*(G0+lambda)+B1[13]*G0*B1[12]+B1[21]*G0*B1[20])*DETJ;
				   Matrix[124] += (B1[21]*G0*B1[3]+B1[5]*lambda*B1[19])*DETJ;
				   Matrix[125] += (B1[13]*G0*B1[3]+B1[5]*lambda*B1[11])*DETJ;
				   Matrix[126] += (B1[5]*B1[3]*(G0+lambda)+B1[13]*G0*B1[11]+B1[21]*G0*B1[19])*DETJ;
				   Matrix[127] += (B1[21]*G0*B1[2]+B1[5]*lambda*B1[18])*DETJ;
				   Matrix[128] += (B1[13]*G0*B1[2]+B1[5]*lambda*B1[10])*DETJ;
				   Matrix[129] += (B1[5]*B1[2]*(G0+lambda)+B1[13]*G0*B1[10]+B1[21]*G0*B1[18])*DETJ;
				   Matrix[130] += (B1[21]*G0*B1[1]+B1[5]*lambda*B1[17])*DETJ;
				   Matrix[131] += (B1[13]*G0*B1[1]+B1[5]*lambda*B1[9])*DETJ;
				   Matrix[132] += (B1[5]*B1[1]*(G0+lambda)+B1[13]*G0*B1[9]+B1[21]*G0*B1[17])*DETJ;
				   Matrix[133] += (B1[21]*G0*B1[0]+B1[5]*lambda*B1[16])*DETJ;
				   Matrix[134] += (B1[13]*G0*B1[0]+B1[5]*lambda*B1[8])*DETJ;
				   Matrix[135] += (B1[5]*B1[0]*(G0+lambda)+B1[13]*G0*B1[8]+B1[21]*G0*B1[16])*DETJ;
				   Matrix[136] += (B1[13]*B1[13]*(G0+lambda)+B1[5]*G0*B1[5]+B1[21]*G0*B1[21])*DETJ;
				   Matrix[137] += (B1[5]*G0*B1[13]+B1[13]*lambda*B1[5])*DETJ;
				   Matrix[138] += (B1[21]*G0*B1[12]+B1[13]*lambda*B1[20])*DETJ;
				   Matrix[139] += (B1[13]*B1[12]*(G0+lambda)+B1[5]*G0*B1[4]+B1[21]*G0*B1[20])*DETJ;
				   Matrix[140] += (B1[5]*G0*B1[12]+B1[13]*lambda*B1[4])*DETJ;
				   Matrix[141] += (B1[21]*G0*B1[11]+B1[13]*lambda*B1[19])*DETJ;
				   Matrix[142] += (B1[13]*B1[11]*(G0+lambda)+B1[5]*G0*B1[3]+B1[21]*G0*B1[19])*DETJ;
				   Matrix[143] += (B1[5]*G0*B1[11]+B1[13]*lambda*B1[3])*DETJ;
				   Matrix[144] += (B1[21]*G0*B1[10]+B1[13]*lambda*B1[18])*DETJ;
				   Matrix[145] += (B1[13]*B1[10]*(G0+lambda)+B1[5]*G0*B1[2]+B1[21]*G0*B1[18])*DETJ;
				   Matrix[146] += (B1[5]*G0*B1[10]+B1[13]*lambda*B1[2])*DETJ;
				   Matrix[147] += (B1[21]*G0*B1[9]+B1[13]*lambda*B1[17])*DETJ;
				   Matrix[148] += (B1[13]*B1[9]*(G0+lambda)+B1[5]*G0*B1[1]+B1[21]*G0*B1[17])*DETJ;
				   Matrix[149] += (B1[5]*G0*B1[9]+B1[13]*lambda*B1[1])*DETJ;
				   Matrix[150] += (B1[21]*G0*B1[8]+B1[13]*lambda*B1[16])*DETJ;
				   Matrix[151] += (B1[13]*B1[8]*(G0+lambda)+B1[5]*G0*B1[0]+B1[21]*G0*B1[16])*DETJ;
				   Matrix[152] += (B1[5]*G0*B1[8]+B1[13]*lambda*B1[0])*DETJ;
				   Matrix[153] += (B1[21]*B1[21]*(G0+lambda)+B1[5]*G0*B1[5]+B1[13]*G0*B1[13])*DETJ;
				   Matrix[154] += (B1[13]*G0*B1[21]+B1[21]*lambda*B1[13])*DETJ;
				   Matrix[155] += (B1[5]*G0*B1[21]+B1[21]*lambda*B1[5])*DETJ;
				   Matrix[156] += (B1[21]*B1[20]*(G0+lambda)+B1[5]*G0*B1[4]+B1[13]*G0*B1[12])*DETJ;
				   Matrix[157] += (B1[13]*G0*B1[20]+B1[21]*lambda*B1[12])*DETJ;
				   Matrix[158] += (B1[5]*G0*B1[20]+B1[21]*lambda*B1[4])*DETJ;
				   Matrix[159] += (B1[21]*B1[19]*(G0+lambda)+B1[5]*G0*B1[3]+B1[13]*G0*B1[11])*DETJ;
				   Matrix[160] += (B1[13]*G0*B1[19]+B1[21]*lambda*B1[11])*DETJ;
				   Matrix[161] += (B1[5]*G0*B1[19]+B1[21]*lambda*B1[3])*DETJ;
				   Matrix[162] += (B1[21]*B1[18]*(G0+lambda)+B1[5]*G0*B1[2]+B1[13]*G0*B1[10])*DETJ;
				   Matrix[163] += (B1[13]*G0*B1[18]+B1[21]*lambda*B1[10])*DETJ;
				   Matrix[164] += (B1[5]*G0*B1[18]+B1[21]*lambda*B1[2])*DETJ;
				   Matrix[165] += (B1[21]*B1[17]*(G0+lambda)+B1[5]*G0*B1[1]+B1[13]*G0*B1[9])*DETJ;
				   Matrix[166] += (B1[13]*G0*B1[17]+B1[21]*lambda*B1[9])*DETJ;
				   Matrix[167] += (B1[5]*G0*B1[17]+B1[21]*lambda*B1[1])*DETJ;
				   Matrix[168] += (B1[21]*B1[16]*(G0+lambda)+B1[5]*G0*B1[0]+B1[13]*G0*B1[8])*DETJ;
				   Matrix[169] += (B1[13]*G0*B1[16]+B1[21]*lambda*B1[8])*DETJ;
				   Matrix[170] += (B1[5]*G0*B1[16]+B1[21]*lambda*B1[0])*DETJ;
				   Matrix[171] += (B1[6]*B1[6]*(G0+lambda)+B1[14]*G0*B1[14]+B1[22]*G0*B1[22])*DETJ;
				   Matrix[172] += (B1[22]*G0*B1[5]+B1[6]*lambda*B1[21])*DETJ;
				   Matrix[173] += (B1[14]*G0*B1[5]+B1[6]*lambda*B1[13])*DETJ;
				   Matrix[174] += (B1[6]*B1[5]*(G0+lambda)+B1[14]*G0*B1[13]+B1[22]*G0*B1[21])*DETJ;
				   Matrix[175] += (B1[22]*G0*B1[4]+B1[6]*lambda*B1[20])*DETJ;
				   Matrix[176] += (B1[14]*G0*B1[4]+B1[6]*lambda*B1[12])*DETJ;
				   Matrix[177] += (B1[6]*B1[4]*(G0+lambda)+B1[14]*G0*B1[12]+B1[22]*G0*B1[20])*DETJ;
				   Matrix[178] += (B1[22]*G0*B1[3]+B1[6]*lambda*B1[19])*DETJ;
				   Matrix[179] += (B1[14]*G0*B1[3]+B1[6]*lambda*B1[11])*DETJ;
				   Matrix[180] += (B1[6]*B1[3]*(G0+lambda)+B1[14]*G0*B1[11]+B1[22]*G0*B1[19])*DETJ;
				   Matrix[181] += (B1[22]*G0*B1[2]+B1[6]*lambda*B1[18])*DETJ;
				   Matrix[182] += (B1[14]*G0*B1[2]+B1[6]*lambda*B1[10])*DETJ;
				   Matrix[183] += (B1[6]*B1[2]*(G0+lambda)+B1[14]*G0*B1[10]+B1[22]*G0*B1[18])*DETJ;
				   Matrix[184] += (B1[22]*G0*B1[1]+B1[6]*lambda*B1[17])*DETJ;
				   Matrix[185] += (B1[14]*G0*B1[1]+B1[6]*lambda*B1[9])*DETJ;
				   Matrix[186] += (B1[6]*B1[1]*(G0+lambda)+B1[14]*G0*B1[9]+B1[22]*G0*B1[17])*DETJ;
				   Matrix[187] += (B1[22]*G0*B1[0]+B1[6]*lambda*B1[16])*DETJ;
				   Matrix[188] += (B1[14]*G0*B1[0]+B1[6]*lambda*B1[8])*DETJ;
				   Matrix[189] += (B1[6]*B1[0]*(G0+lambda)+B1[14]*G0*B1[8]+B1[22]*G0*B1[16])*DETJ;
				   Matrix[190] += (B1[14]*B1[14]*(G0+lambda)+B1[6]*G0*B1[6]+B1[22]*G0*B1[22])*DETJ;
				   Matrix[191] += (B1[6]*G0*B1[14]+B1[14]*lambda*B1[6])*DETJ;
				   Matrix[192] += (B1[22]*G0*B1[13]+B1[14]*lambda*B1[21])*DETJ;
				   Matrix[193] += (B1[14]*B1[13]*(G0+lambda)+B1[6]*G0*B1[5]+B1[22]*G0*B1[21])*DETJ;
				   Matrix[194] += (B1[6]*G0*B1[13]+B1[14]*lambda*B1[5])*DETJ;
				   Matrix[195] += (B1[22]*G0*B1[12]+B1[14]*lambda*B1[20])*DETJ;
				   Matrix[196] += (B1[14]*B1[12]*(G0+lambda)+B1[6]*G0*B1[4]+B1[22]*G0*B1[20])*DETJ;
				   Matrix[197] += (B1[6]*G0*B1[12]+B1[14]*lambda*B1[4])*DETJ;
				   Matrix[198] += (B1[22]*G0*B1[11]+B1[14]*lambda*B1[19])*DETJ;
				   Matrix[199] += (B1[14]*B1[11]*(G0+lambda)+B1[6]*G0*B1[3]+B1[22]*G0*B1[19])*DETJ;
				   Matrix[200] += (B1[6]*G0*B1[11]+B1[14]*lambda*B1[3])*DETJ;
				   Matrix[201] += (B1[22]*G0*B1[10]+B1[14]*lambda*B1[18])*DETJ;
				   Matrix[202] += (B1[14]*B1[10]*(G0+lambda)+B1[6]*G0*B1[2]+B1[22]*G0*B1[18])*DETJ;
				   Matrix[203] += (B1[6]*G0*B1[10]+B1[14]*lambda*B1[2])*DETJ;
				   Matrix[204] += (B1[22]*G0*B1[9]+B1[14]*lambda*B1[17])*DETJ;
				   Matrix[205] += (B1[14]*B1[9]*(G0+lambda)+B1[6]*G0*B1[1]+B1[22]*G0*B1[17])*DETJ;
				   Matrix[206] += (B1[6]*G0*B1[9]+B1[14]*lambda*B1[1])*DETJ;
				   Matrix[207] += (B1[22]*G0*B1[8]+B1[14]*lambda*B1[16])*DETJ;
				   Matrix[208] += (B1[14]*B1[8]*(G0+lambda)+B1[6]*G0*B1[0]+B1[22]*G0*B1[16])*DETJ;
				   Matrix[209] += (B1[6]*G0*B1[8]+B1[14]*lambda*B1[0])*DETJ;
				   Matrix[210] += (B1[22]*B1[22]*(G0+lambda)+B1[6]*G0*B1[6]+B1[14]*G0*B1[14])*DETJ;
				   Matrix[211] += (B1[14]*G0*B1[22]+B1[22]*lambda*B1[14])*DETJ;
				   Matrix[212] += (B1[6]*G0*B1[22]+B1[22]*lambda*B1[6])*DETJ;
				   Matrix[213] += (B1[22]*B1[21]*(G0+lambda)+B1[6]*G0*B1[5]+B1[14]*G0*B1[13])*DETJ;
				   Matrix[214] += (B1[14]*G0*B1[21]+B1[22]*lambda*B1[13])*DETJ;
				   Matrix[215] += (B1[6]*G0*B1[21]+B1[22]*lambda*B1[5])*DETJ;
				   Matrix[216] += (B1[22]*B1[20]*(G0+lambda)+B1[6]*G0*B1[4]+B1[14]*G0*B1[12])*DETJ;
				   Matrix[217] += (B1[14]*G0*B1[20]+B1[22]*lambda*B1[12])*DETJ;
				   Matrix[218] += (B1[6]*G0*B1[20]+B1[22]*lambda*B1[4])*DETJ;
				   Matrix[219] += (B1[22]*B1[19]*(G0+lambda)+B1[6]*G0*B1[3]+B1[14]*G0*B1[11])*DETJ;
				   Matrix[220] += (B1[14]*G0*B1[19]+B1[22]*lambda*B1[11])*DETJ;
				   Matrix[221] += (B1[6]*G0*B1[19]+B1[22]*lambda*B1[3])*DETJ;
				   Matrix[222] += (B1[22]*B1[18]*(G0+lambda)+B1[6]*G0*B1[2]+B1[14]*G0*B1[10])*DETJ;
				   Matrix[223] += (B1[14]*G0*B1[18]+B1[22]*lambda*B1[10])*DETJ;
				   Matrix[224] += (B1[6]*G0*B1[18]+B1[22]*lambda*B1[2])*DETJ;
				   Matrix[225] += (B1[22]*B1[17]*(G0+lambda)+B1[6]*G0*B1[1]+B1[14]*G0*B1[9])*DETJ;
				   Matrix[226] += (B1[14]*G0*B1[17]+B1[22]*lambda*B1[9])*DETJ;
				   Matrix[227] += (B1[6]*G0*B1[17]+B1[22]*lambda*B1[1])*DETJ;
				   Matrix[228] += (B1[22]*B1[16]*(G0+lambda)+B1[6]*G0*B1[0]+B1[14]*G0*B1[8])*DETJ;
				   Matrix[229] += (B1[14]*G0*B1[16]+B1[22]*lambda*B1[8])*DETJ;
				   Matrix[230] += (B1[6]*G0*B1[16]+B1[22]*lambda*B1[0])*DETJ;
				   Matrix[231] += (B1[7]*B1[7]*(G0+lambda)+B1[15]*G0*B1[15]+B1[23]*G0*B1[23])*DETJ;
				   Matrix[232] += (B1[23]*G0*B1[6]+B1[7]*lambda*B1[22])*DETJ;
				   Matrix[233] += (B1[15]*G0*B1[6]+B1[7]*lambda*B1[14])*DETJ;
				   Matrix[234] += (B1[7]*B1[6]*(G0+lambda)+B1[15]*G0*B1[14]+B1[23]*G0*B1[22])*DETJ;
				   Matrix[235] += (B1[23]*G0*B1[5]+B1[7]*lambda*B1[21])*DETJ;
				   Matrix[236] += (B1[15]*G0*B1[5]+B1[7]*lambda*B1[13])*DETJ;
				   Matrix[237] += (B1[7]*B1[5]*(G0+lambda)+B1[15]*G0*B1[13]+B1[23]*G0*B1[21])*DETJ;
				   Matrix[238] += (B1[23]*G0*B1[4]+B1[7]*lambda*B1[20])*DETJ;
				   Matrix[239] += (B1[15]*G0*B1[4]+B1[7]*lambda*B1[12])*DETJ;
				   Matrix[240] += (B1[7]*B1[4]*(G0+lambda)+B1[15]*G0*B1[12]+B1[23]*G0*B1[20])*DETJ;
				   Matrix[241] += (B1[23]*G0*B1[3]+B1[7]*lambda*B1[19])*DETJ;
				   Matrix[242] += (B1[15]*G0*B1[3]+B1[7]*lambda*B1[11])*DETJ;
				   Matrix[243] += (B1[7]*B1[3]*(G0+lambda)+B1[15]*G0*B1[11]+B1[23]*G0*B1[19])*DETJ;
				   Matrix[244] += (B1[23]*G0*B1[2]+B1[7]*lambda*B1[18])*DETJ;
				   Matrix[245] += (B1[15]*G0*B1[2]+B1[7]*lambda*B1[10])*DETJ;
				   Matrix[246] += (B1[7]*B1[2]*(G0+lambda)+B1[15]*G0*B1[10]+B1[23]*G0*B1[18])*DETJ;
				   Matrix[247] += (B1[23]*G0*B1[1]+B1[7]*lambda*B1[17])*DETJ;
				   Matrix[248] += (B1[15]*G0*B1[1]+B1[7]*lambda*B1[9])*DETJ;
				   Matrix[249] += (B1[7]*B1[1]*(G0+lambda)+B1[15]*G0*B1[9]+B1[23]*G0*B1[17])*DETJ;
				   Matrix[250] += (B1[23]*G0*B1[0]+B1[7]*lambda*B1[16])*DETJ;
				   Matrix[251] += (B1[15]*G0*B1[0]+B1[7]*lambda*B1[8])*DETJ;
				   Matrix[252] += (B1[7]*B1[0]*(G0+lambda)+B1[15]*G0*B1[8]+B1[23]*G0*B1[16])*DETJ;
				   Matrix[253] += (B1[15]*B1[15]*(G0+lambda)+B1[7]*G0*B1[7]+B1[23]*G0*B1[23])*DETJ;
				   Matrix[254] += (B1[7]*G0*B1[15]+B1[15]*lambda*B1[7])*DETJ;
				   Matrix[255] += (B1[23]*G0*B1[14]+B1[15]*lambda*B1[22])*DETJ;
				   Matrix[256] += (B1[15]*B1[14]*(G0+lambda)+B1[7]*G0*B1[6]+B1[23]*G0*B1[22])*DETJ;
				   Matrix[257] += (B1[7]*G0*B1[14]+B1[15]*lambda*B1[6])*DETJ;
				   Matrix[258] += (B1[23]*G0*B1[13]+B1[15]*lambda*B1[21])*DETJ;
				   Matrix[259] += (B1[15]*B1[13]*(G0+lambda)+B1[7]*G0*B1[5]+B1[23]*G0*B1[21])*DETJ;
				   Matrix[260] += (B1[7]*G0*B1[13]+B1[15]*lambda*B1[5])*DETJ;
				   Matrix[261] += (B1[23]*G0*B1[12]+B1[15]*lambda*B1[20])*DETJ;
				   Matrix[262] += (B1[15]*B1[12]*(G0+lambda)+B1[7]*G0*B1[4]+B1[23]*G0*B1[20])*DETJ;
				   Matrix[263] += (B1[7]*G0*B1[12]+B1[15]*lambda*B1[4])*DETJ;
				   Matrix[264] += (B1[23]*G0*B1[11]+B1[15]*lambda*B1[19])*DETJ;
				   Matrix[265] += (B1[15]*B1[11]*(G0+lambda)+B1[7]*G0*B1[3]+B1[23]*G0*B1[19])*DETJ;
				   Matrix[266] += (B1[7]*G0*B1[11]+B1[15]*lambda*B1[3])*DETJ;
				   Matrix[267] += (B1[23]*G0*B1[10]+B1[15]*lambda*B1[18])*DETJ;
				   Matrix[268] += (B1[15]*B1[10]*(G0+lambda)+B1[7]*G0*B1[2]+B1[23]*G0*B1[18])*DETJ;
				   Matrix[269] += (B1[7]*G0*B1[10]+B1[15]*lambda*B1[2])*DETJ;
				   Matrix[270] += (B1[23]*G0*B1[9]+B1[15]*lambda*B1[17])*DETJ;
				   Matrix[271] += (B1[15]*B1[9]*(G0+lambda)+B1[7]*G0*B1[1]+B1[23]*G0*B1[17])*DETJ;
				   Matrix[272] += (B1[7]*G0*B1[9]+B1[15]*lambda*B1[1])*DETJ;
				   Matrix[273] += (B1[23]*G0*B1[8]+B1[15]*lambda*B1[16])*DETJ;
				   Matrix[274] += (B1[15]*B1[8]*(G0+lambda)+B1[7]*G0*B1[0]+B1[23]*G0*B1[16])*DETJ;
				   Matrix[275] += (B1[7]*G0*B1[8]+B1[15]*lambda*B1[0])*DETJ;
				   Matrix[276] += (B1[23]*B1[23]*(G0+lambda)+B1[7]*G0*B1[7]+B1[15]*G0*B1[15])*DETJ;
				   Matrix[277] += (B1[15]*G0*B1[23]+B1[23]*lambda*B1[15])*DETJ;
				   Matrix[278] += (B1[7]*G0*B1[23]+B1[23]*lambda*B1[7])*DETJ;
				   Matrix[279] += (B1[23]*B1[22]*(G0+lambda)+B1[7]*G0*B1[6]+B1[15]*G0*B1[14])*DETJ;
				   Matrix[280] += (B1[15]*G0*B1[22]+B1[23]*lambda*B1[14])*DETJ;
				   Matrix[281] += (B1[7]*G0*B1[22]+B1[23]*lambda*B1[6])*DETJ;
				   Matrix[282] += (B1[23]*B1[21]*(G0+lambda)+B1[7]*G0*B1[5]+B1[15]*G0*B1[13])*DETJ;
				   Matrix[283] += (B1[15]*G0*B1[21]+B1[23]*lambda*B1[13])*DETJ;
				   Matrix[284] += (B1[7]*G0*B1[21]+B1[23]*lambda*B1[5])*DETJ;
				   Matrix[285] += (B1[23]*B1[20]*(G0+lambda)+B1[7]*G0*B1[4]+B1[15]*G0*B1[12])*DETJ;
				   Matrix[286] += (B1[15]*G0*B1[20]+B1[23]*lambda*B1[12])*DETJ;
				   Matrix[287] += (B1[7]*G0*B1[20]+B1[23]*lambda*B1[4])*DETJ;
				   Matrix[288] += (B1[23]*B1[19]*(G0+lambda)+B1[7]*G0*B1[3]+B1[15]*G0*B1[11])*DETJ;
				   Matrix[289] += (B1[15]*G0*B1[19]+B1[23]*lambda*B1[11])*DETJ;
				   Matrix[290] += (B1[7]*G0*B1[19]+B1[23]*lambda*B1[3])*DETJ;
				   Matrix[291] += (B1[23]*B1[18]*(G0+lambda)+B1[7]*G0*B1[2]+B1[15]*G0*B1[10])*DETJ;
				   Matrix[292] += (B1[15]*G0*B1[18]+B1[23]*lambda*B1[10])*DETJ;
				   Matrix[293] += (B1[7]*G0*B1[18]+B1[23]*lambda*B1[2])*DETJ;
				   Matrix[294] += (B1[23]*B1[17]*(G0+lambda)+B1[7]*G0*B1[1]+B1[15]*G0*B1[9])*DETJ;
				   Matrix[295] += (B1[15]*G0*B1[17]+B1[23]*lambda*B1[9])*DETJ;
				   Matrix[296] += (B1[7]*G0*B1[17]+B1[23]*lambda*B1[1])*DETJ;
				   Matrix[297] += (B1[23]*B1[16]*(G0+lambda)+B1[7]*G0*B1[0]+B1[15]*G0*B1[8])*DETJ;
				   Matrix[298] += (B1[15]*G0*B1[16]+B1[23]*lambda*B1[8])*DETJ;
				   Matrix[299] += (B1[7]*G0*B1[16]+B1[23]*lambda*B1[0])*DETJ;
			 }
		}
	}
}

void C8H::ElementGravity(double* bodyforce, double Gravity)
{
	clear(bodyforce,24);
}

//	Calculate element stress 
void C8H::ElementStress(double* stress4, double* Displacement)
{	
	C8HMaterial* material = dynamic_cast<C8HMaterial*>(ElementMaterial);
	double E8H = material->E;
	double nu8H = material->posi_ratio;
	double G0 = E8H/(1+nu8H);
	double lambda = nu8H*E8H/((1+nu8H)*(1-2*nu8H));

	double DISP[24];
	for (unsigned int i = 0; i < 24; i++)
	{
		if (LocationMatrix[i])
		{
			DISP[i] = Displacement[LocationMatrix[i] - 1];
		}
		else
		{
			DISP[i] = 0;
		}
	}

	double GPoint[2];
	GPoint[0] = -1/sqrt(3);
	GPoint[1] = 1/sqrt(3);
	for (unsigned int m = 0; m < 2; m++)
	{
		for (unsigned int n = 0; n < 2; n++)
		{
			for (unsigned int o = 0; o < 2; o++)
			{
				double xi = GPoint[m];
				double eta = GPoint[n];
				double zeta = GPoint[o];
				double GN[12];
				GN[0] = 0.125*(1 - eta)*(1 - zeta);
				GN[1] = 0.125*(1 + eta)*(1 - zeta);
				GN[2] = 0.125*(1 - eta)*(1 + zeta);
				GN[3] = 0.125*(1 + eta)*(1 + zeta);
				GN[4] = 0.125*(1 - xi)*(1 - zeta);
				GN[5] = 0.125*(1 + xi)*(1 - zeta);
				GN[6] = 0.125*(1 - xi)*(1 + zeta);
				GN[7] = 0.125*(1 + xi)*(1 + zeta);
				GN[8] = 0.125*(1 - xi)*(1 - eta);
				GN[9] = 0.125*(1 + xi)*(1 - eta);
				GN[10] = 0.125*(1 + xi)*(1 + eta);
				GN[11] = 0.125*(1 - xi)*(1 + eta);

				double J[9];
				J[0] = GN[0]*(nodes[1]->XYZ[0] - nodes[0]->XYZ[0]) + GN[1]*(nodes[2]->XYZ[0] - nodes[3]->XYZ[0]) + GN[2]*(nodes[5]->XYZ[0] - nodes[4]->XYZ[0]) + GN[3]*(nodes[6]->XYZ[0] - nodes[7]->XYZ[0]);
				J[1] = GN[0]*(nodes[1]->XYZ[1] - nodes[0]->XYZ[1]) + GN[1]*(nodes[2]->XYZ[1] - nodes[3]->XYZ[1]) + GN[2]*(nodes[5]->XYZ[1] - nodes[4]->XYZ[1]) + GN[3]*(nodes[6]->XYZ[1] - nodes[7]->XYZ[1]);
				J[2] = GN[0]*(nodes[1]->XYZ[2] - nodes[0]->XYZ[2]) + GN[1]*(nodes[2]->XYZ[2] - nodes[3]->XYZ[2]) + GN[2]*(nodes[5]->XYZ[2] - nodes[4]->XYZ[2]) + GN[3]*(nodes[6]->XYZ[2] - nodes[7]->XYZ[2]);
				J[3] = GN[4]*(nodes[3]->XYZ[0] - nodes[0]->XYZ[0]) + GN[5]*(nodes[2]->XYZ[0] - nodes[1]->XYZ[0]) + GN[6]*(nodes[7]->XYZ[0] - nodes[4]->XYZ[0]) + GN[7]*(nodes[6]->XYZ[0] - nodes[5]->XYZ[0]);
				J[4] = GN[4]*(nodes[3]->XYZ[1] - nodes[0]->XYZ[1]) + GN[5]*(nodes[2]->XYZ[1] - nodes[1]->XYZ[1]) + GN[6]*(nodes[7]->XYZ[1] - nodes[4]->XYZ[1]) + GN[7]*(nodes[6]->XYZ[1] - nodes[5]->XYZ[1]);
				J[5] = GN[4]*(nodes[3]->XYZ[2] - nodes[0]->XYZ[2]) + GN[5]*(nodes[2]->XYZ[2] - nodes[1]->XYZ[2]) + GN[6]*(nodes[7]->XYZ[2] - nodes[4]->XYZ[2]) + GN[7]*(nodes[6]->XYZ[2] - nodes[5]->XYZ[2]);
				J[6] = GN[8]*(nodes[4]->XYZ[0] - nodes[0]->XYZ[0]) + GN[9]*(nodes[5]->XYZ[0] - nodes[1]->XYZ[0]) + GN[10]*(nodes[6]->XYZ[0] - nodes[2]->XYZ[0]) + GN[11]*(nodes[7]->XYZ[0] - nodes[3]->XYZ[0]);
				J[7] = GN[8]*(nodes[4]->XYZ[1] - nodes[0]->XYZ[1]) + GN[9]*(nodes[5]->XYZ[1] - nodes[1]->XYZ[1]) + GN[10]*(nodes[6]->XYZ[1] - nodes[2]->XYZ[1]) + GN[11]*(nodes[7]->XYZ[1] - nodes[3]->XYZ[1]);
				J[8] = GN[8]*(nodes[4]->XYZ[2] - nodes[0]->XYZ[2]) + GN[9]*(nodes[5]->XYZ[2] - nodes[1]->XYZ[2]) + GN[10]*(nodes[6]->XYZ[2] - nodes[2]->XYZ[2]) + GN[11]*(nodes[7]->XYZ[2] - nodes[3]->XYZ[2]);
				double DETJ = J[0]*J[4]*J[8] - J[0]*J[5]*J[7] - J[1]*J[3]*J[8] + J[1]*J[5]*J[6] + J[2]*J[3]*J[7] - J[2]*J[4]*J[6];

				double INVJ[9];
				INVJ[0] = (J[4]*J[8]-J[5]*J[7])/DETJ;
			    INVJ[1] = -(J[1]*J[8]-J[2]*J[7])/DETJ;
			    INVJ[2] = (J[1]*J[5]-J[2]*J[4])/DETJ;
			    INVJ[3] = -(J[3]*J[8]-J[5]*J[6])/DETJ;
			    INVJ[4] = (J[0]*J[8]-J[2]*J[6])/DETJ;
			    INVJ[5] = -(J[0]*J[5]-J[2]*J[3])/DETJ;
			    INVJ[6] = (J[3]*J[7]-J[4]*J[6])/DETJ;
			    INVJ[7] = -(J[0]*J[7]-J[1]*J[6])/DETJ;
			    INVJ[8] = (J[0]*J[4]-J[1]*J[3])/DETJ;

				double B1[24];
				B1[0] = -GN[0]*INVJ[0]-GN[4]*INVJ[1]-GN[8]*INVJ[2];
			    B1[1] = GN[0]*INVJ[0]-GN[5]*INVJ[1]-GN[9]*INVJ[2];
			    B1[2] = GN[1]*INVJ[0]+GN[5]*INVJ[1]-GN[10]*INVJ[2];
			    B1[3] = -GN[1]*INVJ[0]+GN[4]*INVJ[1]-GN[11]*INVJ[2];
				B1[4] = -GN[2]*INVJ[0]-GN[6]*INVJ[1]+GN[8]*INVJ[2];
				B1[5] = GN[2]*INVJ[0]-GN[7]*INVJ[1]+GN[9]*INVJ[2];
				B1[6] = GN[3]*INVJ[0]+GN[7]*INVJ[1]+GN[10]*INVJ[2];
				B1[7] = -GN[3]*INVJ[0]+GN[6]*INVJ[1]+GN[11]*INVJ[2];
			    B1[8] = -GN[0]*INVJ[3]-GN[4]*INVJ[4]-GN[8]*INVJ[5];
			    B1[9] = GN[0]*INVJ[3]-GN[5]*INVJ[4]-GN[9]*INVJ[5];
			    B1[10] = GN[1]*INVJ[3]+GN[5]*INVJ[4]-GN[10]*INVJ[5];
			    B1[11] = -GN[1]*INVJ[3]+GN[4]*INVJ[4]-GN[11]*INVJ[5];
				B1[12] = -GN[2]*INVJ[3]-GN[6]*INVJ[4]+GN[8]*INVJ[5];
				B1[13] = GN[2]*INVJ[3]-GN[7]*INVJ[4]+GN[9]*INVJ[5];
				B1[14] = GN[3]*INVJ[3]+GN[7]*INVJ[4]+GN[10]*INVJ[5];
				B1[15] = -GN[3]*INVJ[3]+GN[6]*INVJ[4]+GN[11]*INVJ[5];
				B1[16] = -GN[0]*INVJ[6]-GN[4]*INVJ[7]-GN[8]*INVJ[8];
			    B1[17] = GN[0]*INVJ[6]-GN[5]*INVJ[7]-GN[9]*INVJ[8];
			    B1[18] = GN[1]*INVJ[6]+GN[5]*INVJ[7]-GN[10]*INVJ[8];
			    B1[19] = -GN[1]*INVJ[6]+GN[4]*INVJ[7]-GN[11]*INVJ[8];
				B1[20] = -GN[2]*INVJ[6]-GN[6]*INVJ[7]+GN[8]*INVJ[8];
				B1[21] = GN[2]*INVJ[6]-GN[7]*INVJ[7]+GN[9]*INVJ[8];
				B1[22] = GN[3]*INVJ[6]+GN[7]*INVJ[7]+GN[10]*INVJ[8];
				B1[23] = -GN[3]*INVJ[6]+GN[6]*INVJ[7]+GN[11]*INVJ[8];

				double EPS[6];
				clear(EPS,6);
				EPS[0] = B1[0]*DISP[0]+B1[1]*DISP[3]+B1[2]*DISP[6]+B1[3]*DISP[9]+B1[4]*DISP[12]+B1[5]*DISP[15]+B1[6]*DISP[18]+B1[7]*DISP[21];
		   	    EPS[1] = B1[8]*DISP[1]+B1[9]*DISP[4]+B1[10]*DISP[7]+B1[11]*DISP[10]+B1[12]*DISP[13]+B1[13]*DISP[16]+B1[14]*DISP[19]+B1[15]*DISP[22];
			    EPS[2] = B1[16]*DISP[2]+B1[17]*DISP[5]+B1[18]*DISP[8]+B1[19]*DISP[11]+B1[20]*DISP[14]+B1[21]*DISP[17]+B1[22]*DISP[20]+B1[23]*DISP[23];
			    EPS[3] = B1[8]*DISP[2]+B1[9]*DISP[5]+B1[16]*DISP[1]+B1[10]*DISP[8]+B1[17]*DISP[4]+B1[11]*DISP[11]+B1[18]*DISP[7]+B1[12]*DISP[14]+B1[19]*DISP[10]+B1[13]*DISP[17]+B1[20]*DISP[13]+B1[14]*DISP[20]+B1[21]*DISP[16]+B1[15]*DISP[23]+B1[22]*DISP[19]+B1[23]*DISP[22];
			    EPS[4] = B1[0]*DISP[2]+B1[1]*DISP[5]+B1[2]*DISP[8]+B1[3]*DISP[11]+B1[4]*DISP[14]+B1[5]*DISP[17]+B1[6]*DISP[20]+B1[7]*DISP[23]+B1[16]*DISP[0]+B1[17]*DISP[3]+B1[18]*DISP[6]+B1[19]*DISP[9]+B1[20]*DISP[12]+B1[21]*DISP[15]+B1[22]*DISP[18]+B1[23]*DISP[21];
			    EPS[5] = B1[0]*DISP[1]+B1[1]*DISP[4]+B1[8]*DISP[0]+B1[2]*DISP[7]+B1[9]*DISP[3]+B1[3]*DISP[10]+B1[4]*DISP[13]+B1[5]*DISP[16]+B1[6]*DISP[19]+B1[7]*DISP[22]+B1[10]*DISP[6]+B1[11]*DISP[9]+B1[12]*DISP[12]+B1[13]*DISP[15]+B1[14]*DISP[18]+B1[15]*DISP[21];

				stress4[0] = EPS[1]*lambda+EPS[2]*lambda+EPS[0]*(G0+lambda);
			    stress4[1] = EPS[0]*lambda+EPS[2]*lambda+EPS[1]*(G0+lambda);
			    stress4[2] = EPS[0]*lambda+EPS[1]*lambda+EPS[2]*(G0+lambda);
			    stress4[3] = EPS[3]*G0;
			    stress4[4] = EPS[4]*G0;
			    stress4[5] = EPS[5]*G0;

			}
		}
	}
}
