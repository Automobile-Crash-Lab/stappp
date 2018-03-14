/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <iomanip>
#include <iostream>
#include <algorithm>

#include "Bandwidth.h"
#include "Domain.h"

using namespace std;

//	Clear an array
template <class type> void clear( type* a, unsigned int N )
{
	for (unsigned int i = 0; i < N; i++)
		a[i] = 0;
}

//	Constructor
CBandwidth::CBandwidth()
{
	NEQ = 0;
	Graph = nullptr;
	re = false;
}

//	Desconstructor
CBandwidth::~CBandwidth()
{
	delete [] Graph;
}

// Two compare functions for "sort"
bool CBandwidth::less_vector(const vector<unsigned int> m1, const vector<unsigned int> m2) 
{
	return m2.size() < m1.size();
}
bool CBandwidth::less_degree(const vector<unsigned int> m1, const vector<unsigned int> m2) 
{
	return m1[0] < m2[0];
}

// Get necessary data from Domain
void CBandwidth::InitialBand()
{
	CDomain* FEMData = CDomain::Instance();

	NEQ = FEMData->GetNEQ();
	NUMEG = FEMData->GetNUMEG();
	EleGrpList = FEMData->GetEleGrpList();
}

// Calculate the Graph structure, store in an array of vector
void CBandwidth::CalculateGraph()
{
	Graph = new vector<unsigned int>[NEQ];

// loop for element group
	for(unsigned int i = 0; i < NUMEG; i++)
	{
// loop for elements in a group
		for(unsigned int j = 0;j < EleGrpList[i].GetNUME();j++)
		{
// loop for all the nodes in an element
			for(unsigned int k = 0;k < EleGrpList[i].GetElement(j).NEN; k++)
			{
// loop for all the degree of freedom of a node
				for(unsigned int l = 0;l < EleGrpList[i].GetElement(j).nodes[k]->NDF;l++)
				{
					unsigned int dof1 = EleGrpList[i].GetElement(j).nodes[k]->bcode[l];
					if(dof1 != 0)
					{	
						for(unsigned int m = 0;m < EleGrpList[i].GetElement(j).NEN; m++)
						{
							for(unsigned int n = 0;n < EleGrpList[i].GetElement(j).nodes[m]->NDF;n++)
							{
								unsigned int dof2 = EleGrpList[i].GetElement(j).nodes[m]->bcode[n];
								if(dof2 != 0)
								{
									bool existt = false;
									for (vector<unsigned int>::size_type ix = 0; ix != Graph[dof1-1].size(); ++ix)
									{
										if(dof2 == Graph[dof1-1][ix])
										{
											existt = true;
											break;
										}
									}
									if (!existt)
										Graph[dof1-1].push_back(dof2);
								}
							}
						}
					}
				}
			}
		}
	}
}

// find the Diameter points and their levels
void CBandwidth::Diameter()
{
	vector<vector<unsigned int> > levelv;
	vector<vector<unsigned int> > levelu;
	// find the node with minimum degree
	v = 0;
	for(unsigned int i = 0; i < NEQ; i++)
	{
		if (Graph[i].size() < Graph[v].size())
			v = i;
	}
	v++;
// Calculate the level rooted at v
	Calculatelevel(v, &levelv);
// Find the two points
	while(1)
	{
// loop for the last layer, compare the depth
		depth = levelv.size();
		unsigned int label = 0;
// Put the nodes in the last layer in order according their degree
		sort(levelv.back().begin(),levelv.back().end());
		unsigned int width = INT_MAX;
		for (vector<unsigned int>::size_type ix = 0; ix != levelv.back().size(); ++ix)
		{
			vector<vector<unsigned int> > temp;
			Calculatelevel(levelv.back()[ix], &temp);
			if(temp.size() > depth)
			{
				label = levelv.back()[ix];
				levelv = temp;
				break;
			}
			unsigned int width_t = 0;
			for(vector<unsigned int>::size_type jx = 0; jx != temp.size(); ++jx)
			{
				if(temp[jx].size() >  width_t)
					width_t = temp[jx].size();
			}
			if(width_t < width)
			{
				levelu = temp;
				u = levelv.back()[ix];
			}
		}
		if(label == 0)
			break;
		v = label;
	}
	CalculatePair(&levelv,&levelu);
}

void CBandwidth::CalculatePair(vector<vector<unsigned int> >* levelv,vector<vector<unsigned int> >* levelu)
{
// label the points with 2 numbers
	coorv = new unsigned int[NEQ];
	cooru = new unsigned int[NEQ];
	for(vector<unsigned int>::size_type ix = 0; ix != levelv->size(); ++ix)
		for(vector<unsigned int>::size_type jx = 0; jx != (*levelv)[ix].size(); ++jx)
			coorv[(*levelv)[ix][jx]-1] = ix + 1;
	for(vector<unsigned int>::size_type ix = 0; ix != levelu->size(); ++ix)
		for(vector<unsigned int>::size_type jx = 0; jx != (*levelu)[ix].size(); ++jx)
			cooru[(*levelu)[ix][jx]-1] = depth - ix;
}

// Calculate level rooted at any point
void CBandwidth::Calculatelevel(unsigned int tag,vector<vector<unsigned int> >* levelstore)
{
// an array to storage the color of the nodes: 0 white, 1 gray, 2 black
	unsigned int* color = new unsigned int[NEQ];
	clear(color,NEQ);
// a list to storage the gray points
	list<unsigned int> level;
	vector<unsigned int> temp;
//calculate the level structure rooted at tag
	level.push_back(tag);
	level.push_back(INT_MAX);
	color[tag-1] = 1;
	while(level.size() != 1)
	{
		if (level.front() != INT_MAX)
		{
			for (auto iter = Graph[level.front() - 1].begin(); iter != Graph[level.front() - 1].end(); iter++)
			{
				if (color[*iter-1] == 0)
				{
					level.push_back(*iter);
					color[*iter-1] = 1;
				}
			}
			temp.push_back(level.front());
			color[level.front() - 1] = 2;
			level.pop_front();
		}
		else
		{
			levelstore->push_back(temp);
			temp.clear();
			level.pop_front();
			level.push_back(INT_MAX);
		}
	}
	levelstore->push_back(temp);
}

void CBandwidth::Layer()
{
// put points (i,i) in Ni
	leveln = new vector<vector<unsigned int> >[depth];
// an array to label whether points has been put into a layer
	unsigned int* color = new unsigned int[NEQ];
	clear(color,NEQ);
	for(unsigned int i = 0; i < NEQ; i++)
	{
		if (coorv[i]==cooru[i])
		{
			vector<unsigned int> temp;
			temp.push_back(i+1);
			temp.push_back(Graph[i].size());
			leveln[coorv[i]-1].push_back(temp);
			color[i] = 2;
		}
	}
// store the rest points in a vector of vector
// a list to storage the gray points
	list<unsigned int> gray;
	vector<vector<unsigned int> > restpoints;
	for(unsigned int i = 0; i < NEQ; i++)
	{
		vector<unsigned int> temp;
		gray.clear();
		if (color[i] == 0)
		{
			gray.push_back(i+1);
			color[i] = 1;
			while (gray.size()!=0)
			{
				for (auto iter = Graph[gray.front() - 1].begin(); iter != Graph[gray.front() - 1].end(); iter++)
				{
					if (color[*iter-1] == 0)
					{
						gray.push_back(*iter);
						color[*iter-1] = 1;
					}
				}
				color[gray.front() - 1] = 2;
				temp.push_back(gray.front());
				gray.pop_front();
			}
			restpoints.push_back(temp);
		}
	}
// sort the subgraph accordint their size
	sort(restpoints.begin(),restpoints.end(),less_vector);
// put the rest points into Ni
	// calculate n vector
	for (vector<unsigned int>::size_type ix = 0; ix != restpoints.size(); ++ix)
	{
		unsigned int* n = new unsigned int[depth];
		for (unsigned int ix = 0; ix != depth; ++ix)
		{
			n[ix] = leveln[ix].size();
//			cout << "n"<< ix << '\t'<<n[ix] << endl;
		}
// calculate h vector and l vector
		unsigned int* h = new unsigned int[depth];
		unsigned int* l = new unsigned int[depth];
		for (unsigned int ix = 0; ix != depth; ++ix)
		{
			h[ix] = n[ix];
			l[ix] = n[ix];
		}
		for (vector<unsigned int>::size_type iy = 0; iy != restpoints[ix].size(); ++iy)
		{
			h[coorv[restpoints[ix][iy]-1]-1]++;
			l[cooru[restpoints[ix][iy]-1]-1]++;
		}
//compare h and l
		unsigned int hmax = *max_element(h,h+depth);
		unsigned int lmax = *max_element(l,l+depth);
		if (ix == 0)
			re = true;
		if (hmax < lmax)
		{
// put all nodes in left
			if (ix == 0)
				re = false;
			for (vector<unsigned int>::size_type iy = 0; iy != restpoints[ix].size(); ++iy)
			{
				vector<unsigned int> temp;
				temp.push_back(restpoints[ix][iy]);
				temp.push_back(Graph[restpoints[ix][iy]-1].size());
				leveln[coorv[restpoints[ix][iy]-1]-1].push_back(temp);
			}
		}
		else if (hmax > lmax)
// put all nodes in right
			for (vector<unsigned int>::size_type iy = 0; iy != restpoints[ix].size(); ++iy)
			{
				vector<unsigned int> temp;
				temp.push_back(restpoints[ix][iy]);
				temp.push_back(Graph[restpoints[ix][iy]-1].size());
				leveln[cooru[restpoints[ix][iy]-1]-1].push_back(temp);
				coorv[restpoints[ix][iy]-1] = cooru[restpoints[ix][iy]-1];
			}
		else
		{
			for (vector<unsigned int>::size_type iy = 0; iy != restpoints[ix].size(); ++iy)
			{
				vector<unsigned int> temp;
				temp.push_back(restpoints[ix][iy]);
				leveln[cooru[restpoints[ix][iy]-1]-1].push_back(temp);
				coorv[restpoints[ix][iy]-1] = cooru[restpoints[ix][iy]-1];
				/*if (n[coorv[restpoints[ix][iy]-1]-1] <= n[cooru[restpoints[ix][iy]-1]-1]++)
				{
					vector<unsigned int> temp;
					temp.push_back(restpoints[ix][iy]);
					temp.push_back(Graph[restpoints[ix][iy]-1].size());
					leveln[coorv[restpoints[ix][iy]-1]-1].push_back(temp);
					//n[coorv[restpoints[ix][iy]-1]-1]++;
				}
				else
				{
					vector<unsigned int> temp;
					temp.push_back(restpoints[ix][iy]);
					temp.push_back(Graph[restpoints[ix][iy]-1].size());
					leveln[cooru[restpoints[ix][iy]-1]-1].push_back(temp);
					//n[cooru[restpoints[ix][iy]-1]-1]++;
					coorv[restpoints[ix][iy]-1] = cooru[restpoints[ix][iy]-1];
				}*/
			}
		}

	}
}

void CBandwidth::Renumber(unsigned int* Order)
{
	unsigned int* order = new unsigned int[depth];
	vector<unsigned int>* listorder = new vector<unsigned int>[depth];
//为每一层排序
	unsigned int* color = new unsigned int[NEQ];
	clear(color,NEQ);
	for(unsigned int i = 0; i < depth ;i++)
	{
		sort(leveln[i].begin(),leveln[i].end(),less_degree);;
	}
//颜色和起始点
	list<unsigned int>* gray = new list<unsigned int>[depth];
	if(Graph[v-1].size() <= Graph[u-1].size())
	{
		gray[0].push_back(v);
		color[v-1] = 1;
		for (unsigned int i = 0; i < depth; i++)
			order[i] = i;
		if(!re)
			re = true;
		else
			re = false;
	}
	else
	{
		gray[0].push_back(u);
		color[v-1] = 1;
		for (unsigned int i = 0; i < depth; i++)
			order[i] = depth - i - 1;
		if(re)
			re = true;
		else
			re = false;
	}
// loop for every layer in leveln
	for (unsigned int ix = 0; ix != depth; ++ix)
	{
		//层内循环
		for(auto jx = leveln[order[ix]].begin(); jx != leveln[order[ix]].end(); jx++)
		{
		//把gray 内的搞定
			while (gray[order[ix]].size()!=0)
			{
				vector<vector<unsigned int> > temp1;
				vector<vector<unsigned int> > temp2;
				for (auto iter = Graph[gray[order[ix]].front() - 1].begin(); iter != Graph[gray[order[ix]].front() - 1].end(); iter++)
				{
					if (color[*iter-1] == 0)
					{
						if (coorv[*iter - 1] - 1 == order[ix])
						{
							color[*iter-1] = 1;
							vector<unsigned int> temp;
							temp.push_back(*iter);
							temp.push_back(Graph[*iter-1].size());
							temp1.push_back(temp);
						}
						else if (ix != depth - 1 && coorv[*iter - 1] - 1 == order[ix + 1])
						{
							color[*iter-1] = 1;
							vector<unsigned int> temp;
							temp.push_back(*iter);
							temp.push_back(Graph[*iter-1].size());
							temp2.push_back(temp);
						}
					}
				}
			//排序
				sort(temp1.begin(),temp1.end(),less_degree);
				sort(temp2.begin(),temp2.end(),less_degree);
				color[gray[order[ix]].front() - 1] = 2;
				listorder[order[ix]].push_back(gray[order[ix]].front());
				gray[order[ix]].pop_front();
				//加入暂存区
				for (auto iter = temp1.begin(); iter != temp1.end(); iter++)
					gray[order[ix]].push_back((*iter)[0]);
				if (ix != depth - 1)
					for (auto iter = temp2.begin(); iter != temp2.end(); iter++)
						gray[order[ix+1]].push_back((*iter)[0]);
			}
		//现在gray内没有点了，看一看能不能接着加
			if (color[(*jx)[0]-1] == 0)
				gray[order[ix]].push_back((*jx)[0]);
		}
	}
	//输出
	/*for (unsigned int i = 0; i < depth; i++)
	{
		for (auto j = listorder[i].begin(); j != listorder[i].end(); j++)
			cout << *j << '\t';
		cout << endl;
	}
	for (unsigned int i = 0; i < depth; i++)
	{
		for (auto j = leveln[i].begin(); j != leveln[i].end(); j++)
			cout << (*j)[0] << '\t';
		cout << endl;
	}*/
	//写入Order
	unsigned int count = 1;
	for (unsigned int i = 0; i < depth; i++)
	{
		for (auto j = listorder[i].begin(); j != listorder[i].end(); j++)
		{
			if (re)
				Order[*j - 1] = 1 + NEQ - count;
			else
				Order[*j - 1] = count;
			count ++;
		}
	}
//	for (int i = 0; i < NEQ; i++)
//	{
//		cout << Order[i]<<endl;
//	}
}
