/***************************************************************/
/*  FEM++ ��A C++ finite element method code for teaching      */
/*     Computational Dynamics Laboratory                       */
/*     School of Aerospace Engineering, Tsinghua University    */
/*                                                             */
/*     http://www.comdyn.cn/                                   */
/***************************************************************/

#pragma once

#include "FEM.h"

using namespace std;

// �˵�Ԫ��
class Bar : public Element
{
public:
	Bar();

	friend FileReader;

	virtual void ElementStiffness(double* Matrix);

	virtual void assembly(double* Matrix);

	virtual void ComputeColumnHeight(unsigned int* ColumnHeight);

	virtual unsigned int LocalMatrixSpace();
};

// �˵Ĳ���
// �����˽��������
class BarMaterial : public Material
{
public:
	double Area;
};