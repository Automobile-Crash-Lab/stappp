/************************************************/
/*              FEMCPP ��ѧ����                 */
/*     �廪��ѧ���캽��ѧԺ���㶯��ѧ������     */
/*     Truss.h  �˵�Ԫͷ�ļ�                    */
/*     ʵ������ά�ռ�˵�Ԫ��Ͳ�����Ķ���     */
/*     ���ߣ�����                               */
/*     ����޸ģ�2017/06/02                     */
/************************************************/
#pragma once
#include "FEM.h"
using namespace std;
// �˵�Ԫ��
class Bar : public Element
{
public:
	Bar();
	friend FileReader;
	virtual void LocalStiffness(double* Matrix);
	virtual void Assembly(double* Matrix);
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