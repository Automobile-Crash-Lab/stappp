/************************************************/
/*              FEMCPP ��ѧ����                 */
/*     �廪��ѧ���캽��ѧԺ���㶯��ѧ������     */
/*     Solver.cpp  �������Դ�ļ�               */
/*     ���ߣ�����                               */
/*     ����޸ģ�2017/06/12                     */
/************************************************/
#include "Solver.h"
#include <math.h>
#include <iostream>
#define Error 1e-24  // �����жϸն����Ƿ�����
using namespace std;
int MIN(int I, int J)
{
	if (I < J) return I;
	else return J;
};
int MAX(int I, int J)
{
	if (I > J) return I;
	else return J;
};
Solver::Solver(FEM* FEMData) : FEMData(FEMData) {};
// LDLT�ֽ�
void LDLTSolver::LDLT()
{
	double* K = FEMData->GetStiffnessMatrix();
	unsigned int* Address = FEMData->GetDiagonalAddress();
	unsigned int N = FEMData->GetFreedom();
	for (int j = 0; j < N; j++)      //����ѭ��
	{
		double* Columnj = &K[Address[j] - 1];
		int ColumnNumberj = Address[j + 1] - Address[j];   //���еķ���Ԫ������
		int Heightj = j - ColumnNumberj + 1;               //��j�е���߷���Ԫ��λ��
		for (int i = Heightj; i <= j; i++)                 //�������еķ���Ԫ��ѭ��
		{
			int ColumnNumberi = Address[i + 1] - Address[i];  //��i�еķ���Ԫ������
			double* Columni = &K[Address[i] - 1];
			int CurPostion = Address[j] + j - i - 1;
			double C = 0;
			int Heighti = i - ColumnNumberi + 1;           //��i�е���߷���Ԫ��λ��
			int Height = MAX(Heighti, Heightj);            //HeightΪi,j������͵ĸ߶�
			for (int M = Height; M < i; M++)
			{
				int AddressI = Address[i] + i - M - 1;
				int AddressJ = Address[j] + j - M - 1;
				C += K[AddressI] * K[AddressJ] * K[Address[M] - 1];
			}
			if (i == j)
			{
				K[CurPostion] = K[CurPostion] - C;
				if (abs(K[CurPostion]) < Error)
				{
					cout << "����װ��" << i + 1 << "�����ɶ�ʱ���ն�������" << endl;
					exit(4);
				}
			}
			else K[CurPostion] = (K[CurPostion] - C) / K[Address[i] - 1];
		}
	}
};
// �ⷽ�̣�����λ��
void LDLTSolver::ComputeDisplacement()
{
	double* Force = FEMData->GetForce();        //������
	double* K = FEMData->GetStiffnessMatrix();  //�Ѿ�����LDLT�ֽ�ĸն���
	double* U = FEMData->GetDisplacement();     //λ��
	unsigned int* Address = FEMData->GetDiagonalAddress();  //�Խ�Ԫ��λ��
	unsigned int Freedom = FEMData->GetFreedom(); //���ɶ���
	// L * V = F , V��Fռ��ͬ����λ��
	for (int i = 0; i < Freedom; i++)
	{
		int Height = Address[i + 1] - Address[i]; //��i�з���Ԫ�ظ���
		int CurPos = Address[i + 1] - 2;
		for (int M = i - Height + 1; M < i; M++)
		{
			Force[i] -= K[CurPos] * Force[M];
			CurPos--;
		}
	}
	// D * S = V,  S V Fռ��ͬ����λ��
	for (int i = 0; i < Freedom; i++)
	{
		Force[i] = Force[i] / K[Address[i] - 1];
	}
	// LT * U = V
	for (int i = Freedom - 1; i >= 0; i--)
	{
		double C = 0;
		for (int M = Freedom - 1; M > i; M--)
		{
			int Height = Address[M + 1] - Address[M];
			if (M - Height + 1 <= i) 
				C += K[Address[M] - 1 + M - i] * Force[M];
		}
		Force[i] = Force[i] - C;
	}
	for (int i = 0; i < Freedom; i++) U[i] = Force[i];
};
void LDLTSolver::Solve()
{ 
	Outputter* Output = Outputter::Instance();
	LDLT();
	for (int i = 0; i < FEMData->GetLoadCaseNumber(); i++)
	{
		FEMData->AssemblyForce(i + 1);
		ComputeDisplacement();
		Output->OutputLoadInfo(i + 1);
		Output->OutputDisplacement();
	}
	return; 
};