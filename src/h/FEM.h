#pragma once
/************************************************/
/*              FEMCPP ��ѧ����                 */
/*     �廪��ѧ���캽��ѧԺ���㶯��ѧ������     */
/*     FEM.h  ����Ԫ��������ͷ�ļ�              */
/*     ���ߣ�����                               */
/*     ����޸ģ�2017/06/18                     */
/************************************************/
#include <string>
#include <fstream>
#include <vector>
#include "Outputter.h"
#include "Solver.h"
#define _DEBUG_
using namespace std;
class FileReader;
class FEM;
// ���Ϻͽ�����
/************************/
// ���Ϻͽ������ʵĻ��ֻ࣬�洢����ģ��һ���������µĲ���ͨ���̳�ʵ�֡�
/************************/
class Material
{
public:
	double E;  //����ģ��
};
// �����
// ��ά����࣬�洢��XYZ����ƽ�����ɶȵ���Ϣ
class Node
{
public:
	double XYZ[3];  // XYZ����
	int Fix[3];     // Լ��������1��Լ����0�����ɡ�
	unsigned int Freedom[3];  // �������������ɶȱ��
	static unsigned int Dimension;  // ��̬���������ɶ���
	Node(double X, double Y, double Z);
};
// ��Ԫ��
// ��Ԫ���࣬�ṩ�˵�Ԫ�Ļ�����Ϣ�ͽӿ�
// ʵ�־��嵥ԪӦ�̳д���
class Element
{
protected:
	unsigned int NodeNumber;    //��Ԫ�Ľڵ���
	Node** NodeList;            //��Ԫ�Ľ��ָ������
	Material* ElementMaterial;  //��Ԫ�Ĳ��������
public:
	Element() :NodeNumber(0), NodeList(NULL), ElementMaterial(NULL) {};
	virtual void LocalStiffness(double* Matrix) = 0;  //���㵥Ԫ�ն��󣬴洢��Matrix��
	virtual void Assembly(double* Matrix) = 0;        //����Ԫ�ն�����װ���ܸգ�Matrix�Ĵ洢��ʽ��Ҫ��LocalStiffness����
	virtual void ComputeColumnHeight(unsigned int* ColumnHeight) = 0; //�����иߣ���skyline�洢������ʹ��
	virtual unsigned int LocalMatrixSpace() = 0;     //���ص�Ԫ�ն�����ռ�Ŀռ��С����Ҫ��Matrix�Ĵ洢��ʽ����
	friend FileReader;
};
// LoadData �����ݴ�����Ϣ
struct LoadData
{
	unsigned NodeNumber;   //����
	unsigned Direction;    //���ķ���
	double Force;          //���Ĵ�С
};
// FileReader �ļ�������
class FileReader
{
private:
	ifstream Input;
public:
	FileReader(string InputFile);
	virtual bool ReadFile(FEM* FEMData);
};
// FEM�࣬�������洢����Ԫ�㷨�õ�������
class FEM
{
private:
	static FEM* _instance;
	string Title;                     //����
	int Mode;						  //���ģʽ
	unsigned int NodeNumber;          //�����
	Node* NodeList;                   //���
	unsigned int ElementGroupNumber;  //��Ԫ����
	unsigned int* ElementNumber;      //ÿ����Ԫ��ĵ�Ԫ��
	Element** ElementGroupList;       //��Ԫ
	unsigned int* MaterialNumber;     //ÿ��������Ĳ�����
	Material** MaterialGroupList;     //����
	LoadData** LoadList;              //�غ�
	unsigned int LoadCaseNumber;      //������
	unsigned int* LoadNumber;         //ÿ���������غ�����

	unsigned int Freedom;             //���ɶ���
	double* StiffnessMatrix;          //�ܸն���
	unsigned int* DiagonalAddress;    //�ն���Խ�Ԫ�صĵ�ַ
	double* Displacement;             //λ��
	double* Force;                    //������
	FEM();
public:
	static FEM* Instance();
	friend FileReader;
	friend Outputter;
	inline double* GetStiffnessMatrix() { return StiffnessMatrix; }
	inline unsigned int* GetDiagonalAddress() { return DiagonalAddress; }
	inline unsigned int GetFreedom() { return Freedom; }
	inline double* GetForce() { return Force; }
	inline double* GetDisplacement() { return Displacement; }
	inline unsigned int GetLoadCaseNumber() { return LoadCaseNumber; }
	void GenerateFreedom();           //�������ɶ���
	void AllocateStiffnessMatrix();   //����ն���Ĵ洢�ռ�
	void AssemblyStiffnessMatrix();   //��װ�ն���
	bool AssemblyForce(unsigned int LoadCase);   //��װ��LoadCase�����������ɹ��򷵻�true
	bool Initial(FileReader* Reader); //���ļ���������ʼ������ȡ������Ϣ
#ifdef _DEBUG_
	void Info();  //debugʱ�õ��������������
#endif
};