/***************************************************************/
/*  FEM++ ��A C++ finite element method code for teaching      */
/*     Computational Dynamics Laboratory                       */
/*     School of Aerospace Engineering, Tsinghua University    */
/*                                                             */
/*     http://www.comdyn.cn/                                   */
/***************************************************************/

#pragma once
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

// ������ļ���ʹ�õ���ģʽ�������ڸ������е���
class Outputter
{
private:
	static Outputter* _instance;

	ofstream OutputFile;      //�����ļ������
	Outputter(string FileName);

public:
	static Outputter* Instance(string FileName = " ");

	void OutputLogo();

	void OutputNodeInfo();

	void OutputLoadInfo(int LoadCase);  //�����LoadCase����������Ϣ

	void OutputDisplacement();          //���λ��
};