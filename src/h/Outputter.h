/************************************************/
/*              FEMCPP ��ѧ����                 */
/*     �廪��ѧ���캽��ѧԺ���㶯��ѧ������     */
/*     Outputter.h  �����ͷ�ļ�                */
/*     ���ߣ�����                               */
/*     ����޸ģ�2017/06/04                     */
/********************************************��****/
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