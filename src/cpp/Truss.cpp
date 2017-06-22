/************************************************/
/*              FEMCPP ��ѧ����                 */
/*     �廪��ѧ���캽��ѧԺ���㶯��ѧ������     */
/*     Truss.cpp  �˵�ԪԴ�ļ�                  */
/*     ���ߣ�����                               */
/*     ����޸ģ�2017/06/02                     */
/************************************************/
#include "Truss.h"
#include <math.h>
#include <iostream>
using namespace std;
Bar::Bar()
{
	NodeNumber = 2;
	NodeList = new Node*[2];
	ElementMaterial = NULL;
}
void Bar::ComputeColumnHeight(unsigned int* ColumnHeight)
{
	vector<int> Freedom;
	//�洢���е����ɶ�
	for (int N = 0; N < 2; N++)
	{
		Node* CurNode = NodeList[N];
		for (int D = 0; D < Node::Dimension; D++)
		{
			if (CurNode->Freedom[D]) Freedom.push_back(CurNode->Freedom[D]);
		}
	}
	//�����и�
	//Ѱ�ҵ�Freedom[j]�е���ߵĵ�Ԫ
	for (int i = 0; i < Freedom.size(); i++)
	{
		int I = Freedom[i];
		for (int j = i + 1; j < Freedom.size(); j++)
		{
			int J = Freedom[j];
			// ʹ��I < J
			if (I >= J)
			{
				int temp = I;
				I = J;
				J = temp;
			}
			int Height = J - I;
			if (ColumnHeight[J] < Height) ColumnHeight[J] = Height;
		}
	}
}
//���ص�Ԫ�ն�����ռ�ռ��С
//���ڸ˵�Ԫ�ĵ�Ԫ�ն���Ϊ��������Matrix�Ĵ�СΪ���������21��Ԫ��
//���صķ�ʽ��ȻΪ���д洢
unsigned int Bar::LocalMatrixSpace() { return 21; }
//���㵥Ԫ�ն���
void Bar::LocalStiffness(double* Matrix)
{
	for (int i = 0; i < 21; i++) Matrix[i] = 0;
	BarMaterial* CurMaterial = (BarMaterial*)ElementMaterial;
	double k = CurMaterial->Area * CurMaterial->E;
	// ����˳�
	int XYZ[3];
	int XYZ2[6];    //����XYZ�Ķ�����ֱ�ΪX^2, Y^2, Z^2, XY, YZ, XZ
	for (int i = 0; i < 3; i++)
	{
		XYZ[i] = NodeList[1]->XYZ[i] - NodeList[0]->XYZ[i];
	}	
	XYZ2[0] = XYZ[0] * XYZ[0];
	XYZ2[1] = XYZ[1] * XYZ[1];
	XYZ2[2] = XYZ[2] * XYZ[2];
	XYZ2[3] = XYZ[0] * XYZ[1];
	XYZ2[4] = XYZ[1] * XYZ[2];
	XYZ2[5] = XYZ[0] * XYZ[2];
	double L2 = XYZ2[0] + XYZ2[1] + XYZ2[2];
	double L = sqrt(L2);
	k = k / L / L2;
	// ����ն���
	Matrix[0] = k*XYZ2[0];
	Matrix[1] = k*XYZ2[1];
	Matrix[2] = k*XYZ2[3];
	Matrix[3] = k*XYZ2[2];
	Matrix[4] = k*XYZ2[4];
	Matrix[5] = k*XYZ2[5];
	Matrix[6] = k*XYZ2[0];
	Matrix[7] = -k*XYZ2[5];
	Matrix[8] = -k*XYZ2[3];
	Matrix[9] = -k*XYZ2[0];
	Matrix[10] = k*XYZ2[1];
	Matrix[11] = k*XYZ2[3];
	Matrix[12] = -k*XYZ2[4];
	Matrix[13] = -k*XYZ2[1];
	Matrix[14] = -k*XYZ2[3];
	Matrix[15] = k*XYZ2[2];
	Matrix[16] = k*XYZ2[4];
	Matrix[17] = k*XYZ2[5];
	Matrix[18] = -k*XYZ2[2];
	Matrix[19] = -k*XYZ2[4];
	Matrix[20] = -k*XYZ2[5];
}
//��װ�ܸ�
void Bar::Assembly(double* Matrix)
{
	//���㵥Ԫ�ն���
	LocalStiffness(Matrix);
	//��װ�ܸն���
	unsigned int Freedom[6];
	Freedom[0] = NodeList[0]->Freedom[0];
	Freedom[1] = NodeList[0]->Freedom[1];
	Freedom[2] = NodeList[0]->Freedom[2];
	Freedom[3] = NodeList[1]->Freedom[0];
	Freedom[4] = NodeList[1]->Freedom[1];
	Freedom[5] = NodeList[1]->Freedom[2];
	//������װ
	int K = 0;   //��ǰԪ���ڵ�Ԫ�ն����еı��
	int KK = 0;  //��ǰԪ�����ܸ��еı��
	int I, J;    //���ɶȱ��
	FEM* FEMData = FEM::Instance();
	unsigned int* DiagonalAddress = FEMData->GetDiagonalAddress();
	double* StiffnessMatrix = FEMData->GetStiffnessMatrix();
	for (int j = 0; j < 6; j++)
	{
		J = Freedom[j];
		if (!J) continue;
		K = (1 + j)*j / 2 + 1;  //���жԽ�Ԫ��λ��
		KK = DiagonalAddress[J - 1];
		for (int i = 0; i <= j; i++)
		{
			I = Freedom[i];
			if (!I) continue;
			StiffnessMatrix[KK + J - I - 1] += Matrix[K + j - i - 1];
		}
	}
	return;
}