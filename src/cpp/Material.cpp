/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl 
			 << "    Expected set : " << mset + 1 << endl
			 << "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> Area >>rho;	// Young's modulus and section area

	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << Area << setw(16) << rho<< endl;
}



bool C3TMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cout << "*** Error *** Material sets must be inputted in order !" << endl 
			 << "   Expected set : " << mset + 1 << endl
			 << "   Provided set : " << nset << endl;

		return false;
	}

	Input >> E >>mu;	// Young's modulus and section area

	return true;
}

//	Write material data to Stream OutputFile
void C3TMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << mu << endl;
}

//	Read material data from stream Input
bool C4QMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl 
			 << "   Expected set : " << mset + 1 << endl
			 << "   Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> posi_rate >>rho >> thickness;	// Young's modulus and section area

	return true;
}

//	Write material data to Stream OutputFile
void C4QMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << posi_rate << setw(16) << thickness << endl;
}


//	Read material data from stream Input
bool CShellMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl
			 << "   Expected set : " << mset + 1 << endl
			 << "   Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> posi_rate >>rho >> thickness;	// Young's modulus and section area

	return true;
}

//	Write material data to Stream OutputFile
void CShellMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << posi_rate << setw(16) << thickness << endl;
}

//	Read material data from stream Input
bool CPlateMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl 
			 << "    Expected set : " << mset + 1 << endl
			 << "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> nu >> thick >> rho;	// Young's modulus and section area and thickness.
	
	return true;
}

//	Read material data from stream Input
bool C8HMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl 
			 << "    Expected set : " << mset + 1 << endl
			 << "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> posi_ratio >>rho;	// Young's modulus and section area

	return true;
}

//	Write material data to Stream
void CPlateMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << nu << setw(16) << thick << endl;
}

void C8HMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << posi_ratio << setw(16) << rho<< endl;
}

bool CBeamMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cout << "*** Error *** Material sets must be inputted in order !" << endl 
			 << "   Expected set : " << mset + 1 << endl
			 << "   Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> v >>rho>> Area >> Iy >> Iz >> Ip;	// Young's modulus and section area and ���Ծ�

	return true;
}

//	Write material data to Stream OutputFile
void CBeamMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << v << setw(16) << Area << setw(16) << Iy << setw(16) << Iz << setw(16) << Ip << endl;
}
