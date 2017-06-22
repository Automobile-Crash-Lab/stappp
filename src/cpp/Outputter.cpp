/************************************************/
/*              FEMCPP ��ѧ����                 */
/*     �廪��ѧ���캽��ѧԺ���㶯��ѧ������     */
/*     Outputter.h  �����ͷ�ļ�                */
/*     ���ߣ�����                               */
/*     ����޸ģ�2017/06/04                     */
/************************************************/
#include "FEM.h"
#include "Outputter.h"
#include <iostream>
#include <iomanip>
#include <time.h>
using namespace std;
// ʱ���������
// ���Ը��������output
void dsptime(const struct tm * ptm, ostream& output)
{
	char *pxq[] = { "��","һ","��","��","��","��","��" };
	output << ptm->tm_year + 1900 << "��" << ptm->tm_mon + 1 << "��" << ptm->tm_mday << "�� ";
	output << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec << " ";
	output << " ����" << pxq[ptm->tm_wday] << endl;
}
Outputter* Outputter::_instance = NULL;
// ���캯��
Outputter::Outputter(string FileName)
{
	OutputFile.open(FileName);
	if (!OutputFile) exit(3);
}
// ��������
Outputter* Outputter::Instance(string FileName)
{
	if (!_instance) _instance = new Outputter(FileName);
	return _instance;
}
// ��ӡ�ļ�ͷ
void Outputter::OutputLogo()
{
	cout << "***********************************************************" << endl;
	cout << "xxxxxx  xxxxxx  xxx       xx      xxxxx  xxxxxx   xxxxxx" << endl;
	cout << "xx      xx      xxxx     xxx    xxx      xx   xx  xx   xx" << endl;
	cout << "xx      xx      xxxx     x x    xx       xx    xx xx    xx" << endl;
	cout << "xx      xx      xx x    xx x   xx        xx    xx xx    xx" << endl;
	cout << "xx      xx      xx xx   xx x   xx        xx    xx xx    xx" << endl;
	cout << "xxxxxx  xxxxxx  xx xx   x  x   xx        xx   xx  xx   xx" << endl;
	cout << "xx      xx      xx  xx xx  x   xx        xxxxxx   xxxxxx" << endl;
	cout << "xx      xx      xx  xx x   x   xx        xx       xx" << endl;
	cout << "xx      xx      xx   xxx   x   xx        xx       xx" << endl;
	cout << "xx      xx      xx   xxx   x    xxx   x  xx       xx" << endl;
	cout << "xx      xxxxxxx xx    x    x      xxxxx  xx       xx" << endl;
	cout << "***********************************************************" << endl;
	cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
	OutputFile << "***********************************************************" << endl;
	OutputFile << "xxxxxx  xxxxxx  xxx       xx      xxxxx  xxxxxx   xxxxxx" << endl;
	OutputFile << "xx      xx      xxxx     xxx    xxx      xx   xx  xx   xx" << endl;
	OutputFile << "xx      xx      xxxx     x x    xx       xx    xx xx    xx" << endl;
	OutputFile << "xx      xx      xx x    xx x   xx        xx    xx xx    xx" << endl;
	OutputFile << "xx      xx      xx xx   xx x   xx        xx    xx xx    xx" << endl;
	OutputFile << "xxxxxx  xxxxxx  xx xx   x  x   xx        xx   xx  xx   xx" << endl;
	OutputFile << "xx      xx      xx  xx xx  x   xx        xxxxxx   xxxxxx" << endl;
	OutputFile << "xx      xx      xx  xx x   x   xx        xx       xx" << endl;
	OutputFile << "xx      xx      xx   xxx   x   xx        xx       xx" << endl;
	OutputFile << "xx      xx      xx   xxx   x    xxx   x  xx       xx" << endl;
	OutputFile << "xx      xxxxxxx xx    x    x      xxxxx  xx       xx" << endl;
	OutputFile << "***********************************************************" << endl;
	OutputFile << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
	FEM* FEMData = FEM::Instance();
	cout << "TITLE : " << FEMData->Title << endl;
	OutputFile << "TITLE : " << FEMData->Title << endl;
	time_t nowtime;
	struct tm *local = new struct tm;
	nowtime = time(NULL);
	localtime_s(local, &nowtime);
	dsptime(local, cout);
	dsptime(local, OutputFile);
}
// ��ӡ�����Ϣ
void Outputter::OutputNodeInfo()
{
	FEM* FEMData = FEM::Instance();
	Node* NodeList = FEMData->NodeList;
	int Page = 30;      // ÿҳ����
	cout << setiosflags(ios::scientific);
	OutputFile << setiosflags(ios::scientific);
	cout << "*********************  N O D E **************************" << endl;
	OutputFile << "*********************  N O D E **************************" << endl;
	for (int i = 0; i < FEMData->NodeNumber; i++)
	{
		if (i % Page == 0)
		{
			cout << "No. ......... X ........... Y ........... Z ..........." << endl;
			OutputFile << "No. ......... X ........... Y ........... Z ..........." << endl;
		}
		cout << setw(12) << i + 1 << setw(14) << NodeList[i].XYZ[0] << setw(14) << NodeList[i].XYZ[1] << setw(14) << NodeList[i].XYZ[2] << endl;
		OutputFile << setw(12) << i + 1 << setw(14) << NodeList[i].XYZ[0] << setw(14) << NodeList[i].XYZ[1] << setw(14) << NodeList[i].XYZ[2] << endl;
	}
}
// ��ӡ�غ���Ϣ
void Outputter::OutputLoadInfo(int LoadCase)
{
	FEM* FEMData = FEM::Instance();
	unsigned int* LoadNumber = FEMData->LoadNumber;
	if (LoadCase > FEMData->LoadCaseNumber) return;
	LoadData* Load = FEMData->LoadList[LoadCase - 1];
	cout << setiosflags(ios::scientific);
	OutputFile << setiosflags(ios::scientific);
	cout << "****************** L O A D C A S E " << LoadCase <<  " ********************" << endl;
	OutputFile << "****************** L O A D C A S E " << LoadCase <<  " ********************" << endl;
	cout << "No. ...... NodeNumber .. DIR ....... ....Force .........." << endl;
	OutputFile << "No. ...... NodeNumber .. DIR ....... ....Force .........." << endl;
	for (int i = 0; i < LoadNumber[LoadCase - 1]; i++)
	{
		cout << setw(10) << i + 1 << setw(14) << Load[i].NodeNumber << setw(12) << Load[i].Direction << setw(17) << Load[i].Force << endl;
		OutputFile << setw(10) << i + 1 << setw(14) << Load[i].NodeNumber << setw(12) << Load[i].Direction << setw(17) << Load[i].Force << endl;
	}
}
// ���λ��
void Outputter::OutputDisplacement()
{
	FEM* FEMData = FEM::Instance();
	int Page = 30;
	Node* NodeList = FEMData->NodeList;
	double* Displacement = FEMData->Displacement;
	cout << setiosflags(ios::scientific);
	OutputFile << setiosflags(ios::scientific);
	cout << "************* D I S P L A C E l M E N T *****************" << endl;
	OutputFile << "************* D I S P L A C E l M E N T *****************" << endl;
	for (int i = 0; i < FEMData->NodeNumber; i++)
	{
		if (i % Page == 0)
		{
			cout << "No. ......... X ........... Y ........... Z ..........." << endl;
			OutputFile << "No. ......... X ........... Y ........... Z ..........." << endl;
		}
		cout << setw(12) << i + 1;
		OutputFile << setw(12) << i + 1;
		for (int j = 0; j < 3; j++)
		{
			if (NodeList[i].Freedom[j] == 0)
			{
				cout << setw(14) << 0.0;
				OutputFile << setw(14) << 0.0;
			}
			else
			{
				cout << setw(14) << Displacement[NodeList[i].Freedom[j] - 1];
				OutputFile << setw(14) << Displacement[NodeList[i].Freedom[j] - 1];
			}
		}
		cout << endl;
		OutputFile << endl;
	}
}