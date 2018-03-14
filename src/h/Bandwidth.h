/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include <string>
#include <fstream>
#include <vector>
#include <list> 

#include "Node.h"
#include "ElementGroup.h"

using namespace std;

template <class type> void clear( type* a, unsigned int N );	// Clear an array


class CBandwidth
{
private:

//! store the graph as an array of vector
	vector<unsigned int>* Graph;

//!	Total number of equations in the system, fetch from Domain::FEMData
	unsigned int NEQ;

//!	Total number of Element Group in the system, fetch from Domain::FEMData
	unsigned int NUMEG;

//!	Element group list, fetch from Domain::FEMData
	 CElementGroup* EleGrpList;

//! Two components of the pair numbers
	unsigned int* coorv;
	unsigned int* cooru;

//! The depth of the final level
	unsigned int depth;

//! The final level;
	vector< vector<unsigned int> >* leveln;

//! Two start points of two levels
	unsigned int v;
	unsigned int u;

//! A bool value to determine whether reverse the order;
	bool re;

public:

//!	Constructor
	CBandwidth();

//!	Desconstructor
	~CBandwidth();

	void InitialBand();

//! Calculate the graph
	void CalculateGraph();

//! Calculate diameter points
	void Diameter();

//! Renumber
	void Renumber(unsigned int* Order);

//! Calculate the pair number of nodes
	void CalculatePair(vector<vector<unsigned int> >* levelv,vector<vector<unsigned int> >* levelu);

//! Put nodes into layers
	void Layer();

//! Calculate levels
	void Calculatelevel(unsigned int tag,vector<vector<unsigned int> >* levelstore);

//! Two compare functions for "sort"
	static bool less_vector(const vector<unsigned int> m1, const vector<unsigned int> m2);
	static bool less_degree(const vector<unsigned int> m1, const vector<unsigned int> m2);
};
