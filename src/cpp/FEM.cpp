/************************************************/
/*              FEMCPP ��ѧ����                 */
/*     �廪��ѧ���캽��ѧԺ���㶯��ѧ������     */
/*     FEM.h  ����Ԫ��������Դ�ļ�              */
/*     ���ߣ�����                               */
/*     ����޸ģ�2017/05/27                     */
/************************************************/
#include "FEM.h"
#include "Truss.h"
#include <iomanip>
#include <iostream>
using namespace std;
/***Node***/
// Node���͵����ɶ���Ϊ3
// ��չת�����ɶ���Ҫ���Ĵ���ֵ
unsigned int Node::Dimension = 3;
Node::Node(double X = 0, double Y = 0, double Z = 0) 
{
	XYZ[0] = X;
	XYZ[1] = Y;
	XYZ[2] = Z;
	Fix[0] = 0;
	Fix[1] = 0;
	Fix[2] = 0;
	Freedom[0] = 0;
	Freedom[1] = 0;
	Freedom[2] = 0;
};
/***FEM***/
// Instance��������
FEM* FEM::_instance = NULL;
FEM::FEM() :Title(""), Mode(0), NodeNumber(0), NodeList(NULL), ElementGroupNumber(0),
ElementGroupList(NULL), ElementNumber(NULL),  MaterialNumber(NULL),
MaterialGroupList(NULL), LoadCaseNumber(0), LoadNumber(NULL), LoadList(NULL), Freedom(0), StiffnessMatrix(NULL), Displacement(NULL),
Force(NULL) {};
FEM* FEM::Instance()
{
	if (!_instance) _instance = new FEM;
	return _instance;
}
// Initial���ļ���ʼ��
// ����Reader�Ķ�ȡ���ݺ���
// ��������ԭ���ǿ���ͨ���̳еķ�ʽ��д�ļ���ȡ�������Ӷ����䲻ͬ���ļ���ʽ
bool FEM::Initial(FileReader* Reader)
{
	return Reader->ReadFile(this);
}
/***FileReader***/
FileReader::FileReader(string Filename)
{
	Input.open(Filename);
	if (!Input) exit(3);
}
// ��ȡ�����ļ�
bool FileReader::ReadFile(FEM* FEMData)
{
	int N, LL;
	Input >> FEMData->Title;
	// ���������
	Input >> FEMData->NodeNumber >> FEMData->ElementGroupNumber >> FEMData->LoadCaseNumber
		>> FEMData->Mode;
	// ������
	FEMData->NodeList = new Node[FEMData->NodeNumber];
	Node* NodeList = FEMData->NodeList;
	for (int i = 0; i < FEMData->NodeNumber; i++)
	{
		Input >> N;
		if (N != i + 1) return false;
		Input >> NodeList[i].Fix[0] >> NodeList[i].Fix[1] >> NodeList[i].Fix[2]
			>> NodeList[i].XYZ[0] >> NodeList[i].XYZ[1] >> NodeList[i].XYZ[2];
	}
	// ���빤��
	FEMData->LoadList = new LoadData*[FEMData->LoadCaseNumber];
	FEMData->LoadNumber = new unsigned int[FEMData->LoadCaseNumber];
	LoadData** LoadList = FEMData->LoadList;
	unsigned int* LoadNumber = FEMData->LoadNumber;
	for (int i = 0; i < FEMData->LoadCaseNumber; i++)
	{
		Input >> LL >> N;
		if (LL != i + 1) return false;
		LoadNumber[i] = N;
		LoadList[i] = new LoadData[LoadNumber[i]];
		for (int j = 0; j < LoadNumber[i]; j++)
		{
			Input >>LoadList[i][j].NodeNumber >> LoadList[i][j].Direction >> LoadList[i][j].Force;
		}
	}
	// ���뵥Ԫ������
	FEMData->ElementNumber = new unsigned int[FEMData->ElementGroupNumber];
	unsigned int* ElementNumber = FEMData->ElementNumber;
	FEMData->ElementGroupList = new Element*[FEMData->ElementGroupNumber];
	FEMData->MaterialNumber = new unsigned int[FEMData->ElementGroupNumber];
	unsigned int* MaterialNumber = FEMData->MaterialNumber;
	FEMData->MaterialGroupList = new Material*[FEMData->ElementGroupNumber];
	unsigned int ElementType;
	for (int i = 0; i < FEMData->ElementGroupNumber; i++)
	{
		Input >> ElementType >> ElementNumber[i] >> MaterialNumber[i];
		// ���ݲ�ͬ�ĵ�Ԫ���Ͷ�ȡ���� 
		// �˴����Գ��Ը�д�ø������
		switch (ElementType)
		{
		case 1:
		{
			BarMaterial* MaterialGroup = new BarMaterial[MaterialNumber[i]];
			FEMData->MaterialGroupList[i] = MaterialGroup;
			for (int j = 0; j < MaterialNumber[i]; j++)
			{
				Input >> N;
				if (N != j + 1) return false;
				Input >> MaterialGroup[j].E >> MaterialGroup[j].Area;
			}
			Bar* ElementList = new Bar[ElementNumber[i]];
			FEMData->ElementGroupList[i] = ElementList;
			for (int j = 0; j < ElementNumber[i]; j++)
			{
				Input >> N;
				if (N != j + 1) return false;
				int MNumber;
				int N1, N2;
				Input >> N1 >> N2 >> MNumber;
				ElementList[j].ElementMaterial = &MaterialGroup[MNumber - 1];
				ElementList[j].NodeList[0] = &NodeList[N1 - 1];
				ElementList[j].NodeList[1] = &NodeList[N2 - 1];
			}
			break;
		}
		default:
			return false;
		}
	}
	return true;
}
/***FEM***/
// �������ɶ���
// �㷨�����˴�Լ������Freedom��0
//       �����ɣ��������ɶ� + 1
void FEM::GenerateFreedom()
{
	int CurFreedom = 0;
	for (int i = 0; i < NodeNumber; i++)
	{
		for (int j = 0; j < Node::Dimension; j++)
		{
			if (NodeList[i].Fix[j]) NodeList[i].Freedom[j] = 0;
			else
			{
				CurFreedom++;
				NodeList[i].Freedom[j] = CurFreedom;
			}
		}
	}
	Freedom = CurFreedom;
}
// ����ն���Ĵ洢�ռ䣬ͬʱ����DiagonalAddress����
void FEM::AllocateStiffnessMatrix()
{
	Displacement = new double[Freedom];
	Force = new double[Freedom];
	DiagonalAddress = new unsigned int[Freedom + 1];
	// ����ÿ����Ԫ�Ĳ�ͬ���ɶ�I < J, ˵���ڵ�(I,J)λ����Ԫ�أ���ô��J�е��и�����ΪJ - I
	for (int i = 0; i < Freedom + 1; i++) DiagonalAddress[i] = 0;
	for (int EG = 0; EG < ElementGroupNumber; EG++)
	{
		for (int EN = 0; EN < ElementNumber[EG]; EN++)
		{
			ElementGroupList[EG][EN].ComputeColumnHeight(DiagonalAddress);
		}
	}
	// ���ݼ���õ����и�DiagonalAddress������ռ䣬����Խ�Ԫ��ַ
	DiagonalAddress[0] = 1;
	for (int C = 1; C <= Freedom; C++)
	{
		DiagonalAddress[C] = DiagonalAddress[C - 1] + DiagonalAddress[C] + 1;
	}
	StiffnessMatrix = new double[DiagonalAddress[Freedom] - 1];
	for (int i = 0; i < DiagonalAddress[Freedom] - 1; i++) StiffnessMatrix[i] = 0;
}
//��װ�ն���
//��ÿ����Ԫ�飨ͬ�ֵ�Ԫ������һ��Matrix�ռ�
//����ÿ����Ԫ��Assembly����
//��װ��ʽ�ڵ�Ԫ��ʵ��
void FEM::AssemblyStiffnessMatrix()
{
	for (int EG = 0; EG < ElementGroupNumber; EG++)
	{
		unsigned int Space = ElementGroupList[EG][0].LocalMatrixSpace();
		double* Matrix = new double[Space];
		for (int E = 0; E < ElementNumber[EG]; E++)
		{
			ElementGroupList[EG][E].Assembly(Matrix);
		}
		delete [] Matrix;
	}
}
//��װ������
bool FEM::AssemblyForce(unsigned int LoadCase)
{
	if (LoadCase > LoadCaseNumber) return false;
	LoadData* Load = LoadList[LoadCase - 1];
	//�Ƚ�����Ϊ��
	for (int i = 0; i < Freedom; i++) Force[i] = 0;
	//��װ������
	for (int FN = 0; FN < LoadNumber[LoadCase - 1]; FN++)
	{
		int FreedomDegree = NodeList[Load[FN].NodeNumber - 1].Freedom[Load[FN].Direction - 1];
		Force[FreedomDegree - 1] += Load[FN].Force;
	}
	return true;
}
#ifdef _DEBUG_
//Info
//debugʱ�õ������������
//��ǰ��������ն����һά�洢��ʽ���Խ�Ԫ�ص�λ���Լ���ά�ĸն���
void FEM::Info()
{
	cout << "StiffnessMatrix : " << endl;
	for (int i = 0; i < DiagonalAddress[Freedom] - 1; i++) cout << StiffnessMatrix[i] << " ";
	cout << endl;
	cout << "Address : " << endl;
	for (int i = 0; i < Freedom + 1; i++) cout << DiagonalAddress[i] << " ";
	cout << endl;
	cout << "Matrix : " << endl;
	for (int I = 0; I < Freedom; I++)
	{
		for (int J = 0; J < Freedom; J++)
		{
			int i = I;
			int j = J;
			if (i > j)
			{
				int temp = i;
				i = j;
				j = temp;
			}
			cout << setiosflags(ios::scientific);
			int H = DiagonalAddress[j + 1] - DiagonalAddress[j];
			if (j - i - H >= 0) cout << setw(15) << 0.0;
			else cout << setw(15) << StiffnessMatrix[DiagonalAddress[j] + j - i - 1] << "  ";
		}
		cout << endl;
	}
	cout << endl;
	cout << "U : " << endl;
	for (int I = 0; I < Freedom; I++)
	{
		cout << Displacement[I] << " ";
	}
	cout << endl;
}
#endif