/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.0, October 14, 2017                                         */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once
#include "Outputter.h"
#include <stddef.h>
#include <iostream>
#include <fstream>

using namespace std;

//!	Material base class which only define one data member
/*!	All type of material classes should be derived from this base class */
class CMaterial
{
public:

	unsigned int nset;	//!< Number of set
	
	double E;  //!< Young's modulus

	double rho;
public:

	//! Virtual deconstructor
    virtual ~CMaterial() {};

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset) = 0;

//!	Write material data to Stream
    virtual void Write(COutputter& output, unsigned int mset) = 0;

};

//!	Material class for bar element
class CBarMaterial : public CMaterial
{
public:

	double Area;	//!< Sectional area of a bar element

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream 
	virtual void Write(COutputter& output, unsigned int mset);
};

//!	Material class for 3T element
class C3TMaterial : public CMaterial
{
public:

	double mu;

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream OutputFile
	virtual void Write(COutputter& output, unsigned int mset);
};



//!	Material class for 4Q element
class C4QMaterial : public CMaterial
{
public:

	double posi_rate;	//!< Sectional area of a 4Q element

	double thickness;

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream OutputFile
	virtual void Write(COutputter& output, unsigned int mset);
};

//!	Material class for 4Q element
class CShellMaterial : public CMaterial
{
public:

	double posi_rate;	//!< Sectional area of a 4Q element

	double thickness;
//!	Material class for Plate element
public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);

};


class CPlateMaterial : public CMaterial
{
public:

	double nu;	//!< Possion raito

	double thick;   //!< Thickness of the plate.

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};


class C8HMaterial : public CMaterial
{
public:

	double posi_ratio;	//!< Sectional area of a 8H element

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream OutputFile
	virtual void Write(COutputter& output, unsigned int mset);
};

class CBeamMaterial : public CMaterial
{
public:

	double Iy;	//���Ծ�
	double Iz;	//���Ծ�
	double Ip;	//x������Ծ�
	double v;	//Poisson Ratio
	double Area;	//!< Sectional area of a bar element

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream OutputFile
	virtual void Write(COutputter& output, unsigned int mset);
};
