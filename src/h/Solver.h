/************************************************/
/*              FEMCPP ��ѧ����                 */
/*     �廪��ѧ���캽��ѧԺ���㶯��ѧ������     */
/*     Solver.h  �������ͷ�ļ�                 */
/*     ���ߣ�����                               */
/*     ����޸ģ�2017/06/05                     */
/************************************************/
#pragma once
#include "FEM.h"
using namespace std;
// ������࣬�ṩ�������Ψһ�ӿ�Solve
// ʵ���µ��������Ҫ�̳д���
// ��Ҫ��FEM������ݴ洢��ʽ��ƥ��
class FEM;
class Solver
{
protected:
	FEM* FEMData;
public:
	Solver(FEM* FEMData);
	virtual void Solve() = 0;
};
// LDLT�ֽ������
// ��skyline�洢��ʽƥ��
class LDLTSolver : public Solver
{
public:
	LDLTSolver(FEM* FEMData) :Solver(FEMData) {};
	void LDLT();                // ���ܸն���ִ��LDLT�ֽ⣬�ֽ��洢��ԭλ�� 
	void ComputeDisplacement(); // ʹ�õ�ǰ������������λ�� 
	virtual void Solve();       // ��FEM����
};