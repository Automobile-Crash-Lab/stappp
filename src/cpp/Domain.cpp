/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Material.h"
#include "qsort.h"
#include "Solver_Sparse.h"
#include <iomanip>
#include <iostream>

 template<class T_>
 void debugout(T_* a,int n){
     for(int i=0;i<n;i++){
         cout<< a[i]<<endl;
     }
 }

 void debugoutstiff(double*a){
     for(int i=0;i<12;i++){
         for(int j=0;j<12;j++){
             cout<< setiosflags(ios::scientific) <<setprecision(3) << setw(12)<<((i>j)?a[i*(i+3)/2-j]:a[j*(j+3)/2-i]);
         }
         cout<< endl;
     }
 }

using namespace std;

//	Clear an array
template <class type> void clear( type* a, unsigned int N )
{
	for (unsigned int i = 0; i < N; i++)
		a[i] = 0;
}

CDomain* CDomain::_instance = nullptr;

//	Constructor
CDomain::CDomain()
{
	Title[0] = '0';
	MODEX = 0; //求解模式，1求解，0检查数据

	NUMNP = 0;//节点总数
	NodeList = nullptr;
	
	NUMEG = 0; //单元组总数，每个里面都只能有相同形状的单元，但可以有不同材料、长度。
	EleGrpList = nullptr;
	
	NLCASE = 0; //载荷工况数
	NLOAD = nullptr; //本工况中集中载荷的个数
	LoadCases = nullptr; //工况数编号
	
	NEQ = 0;//总方程数
	NWK = 0;//总刚度阵中元素总数
	MK = 0;//最大列高

	Force = nullptr;
	StiffnessMatrix = nullptr;
	solver=nullptr;
	Gravity=-9.8;
}

//	Desconstructor
CDomain::~CDomain()
{
	delete [] NodeList;

	delete [] EleGrpList;

	delete [] NLOAD;
	delete [] LoadCases;

	delete [] Force;
	delete [] StiffnessMatrix;
}

//	Return pointer to the instance of the Domain class
CDomain* CDomain::Instance()
{
	if (!_instance) 
		_instance = new CDomain();
	
	return _instance;
}

//	Read domain data from the input data file
bool CDomain::ReadData(string FileName, string OutFile)
{
	Input.open(FileName);

	if (!Input) 
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}

	COutputter* Output = COutputter::Instance(OutFile);

//	Read the heading line
	Input.getline(Title, 256);
	Output->OutputHeading();

//	Read the control line
	Input >> NUMNP >> NUMEG >> NLCASE >> MODEX;

//	Read nodal point data
    if (ReadNodalPoints()){
        Output->OutputNodeInfo();
    }
    else
        return false;

//	Update equation number
	CalculateEquationNumber();
	Output->OutputEquationNumber();

//	Read load data
    if (ReadLoadCases()){
        Output->OutputLoadInfo();
    }
    else
        return false;

//	Read element data
    if (ReadElements()){
        Output->OutputElementInfo();
    }
    else
        return false;

	return true;
}

//	Read nodal point data
bool CDomain::ReadNodalPoints()
{

//	Read nodal point data lines
	NodeList = new CNode[NUMNP];

//	Loop over for all nodal points
	for (unsigned int np = 0; np < NUMNP; np++)
		if (!NodeList[np].Read(Input, np))
			return false;

	return true;
}

//	Calculate global equation numbers corresponding to every degree of freedom of each node
void CDomain::CalculateEquationNumber()
{
	NEQ = 0;
	for (unsigned int np = 0; np < NUMNP; np++)	// Loop over for all node
	{
		for (unsigned int dof = 0; dof < 6; dof++)	// Loop over for DOFs of node np
		{
			if (NodeList[np].bcode[dof]) 
				NodeList[np].bcode[dof] = 0;
			else
			{
				NEQ++;
				NodeList[np].bcode[dof] = NEQ;
			}
		}
	}
}

//	Calculate global equation numbers corresponding to every degree of freedom of each node
void CDomain::RefreshEquationNumber()
{
	for (unsigned int np = 0; np < NUMNP; np++)	// Loop over for all node
	{
		for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
		{
			if (NodeList[np].bcode[dof] != 0)
			{
				NodeList[np].bcode[dof] = Order[NodeList[np].bcode[dof]-1];
			}
		}
	}
	delete[] Order;
}

//	Read load case data
bool CDomain::ReadLoadCases()
{
//	Read load data lines
	LoadCases = new CLoadCaseData[NLCASE];	// List all load cases

//	Loop over for all load cases
	for (unsigned int lcase = 0; lcase < NLCASE; lcase++)
		if (!LoadCases[lcase].Read(Input, lcase))
			return false;

	return true;
}

// Read element data
bool CDomain::ReadElements()
{
    EleGrpList = new CElementGroup[NUMEG];

//	Loop over for all element group
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
        if (!EleGrpList[EleGrp].Read(Input))
            return false;
    
    return true;
}

//	Calculate column heights
void CDomain::CalculateColumnHeights()
{    
	unsigned int* ColumnHeights = StiffnessMatrix->GetColumnHeights();

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
    {
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();
        
		for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
			ElementGrp.GetElement(Ele).CalculateColumnHeight(ColumnHeights);
    }

//	Maximum half bandwidth ( = max(ColumnHeights) + 1 )
	MK = ColumnHeights[0];

	for (unsigned int i=1; i<NEQ; i++)
		if (MK < ColumnHeights[i])
			MK = ColumnHeights[i];

	MK = MK + 1;
    
#ifdef _DEBUG_
	COutputter* Output = COutputter::Instance();
	//Output->PrintColumnHeights();
#endif

}

//	Calculate address of diagonal elements in banded matrix
//	Caution: Address is numbered from 1 !
void CDomain::CalculateDiagnoalAddress()
{
    unsigned int* ColumnHeights = StiffnessMatrix->GetColumnHeights();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();
//	clear(DiagonalAddress, NEQ + 1);	// Set all elements to zero

//	Calculate the address of diagonal elements
//	M(0) = 1;  M(i+1) = M(i) + H(i) + 1 (i = 0:NEQ)
	DiagonalAddress[0] = 1;
	for (unsigned int col = 1; col <= NEQ; col++)
		DiagonalAddress[col] = DiagonalAddress[col - 1] + ColumnHeights[col-1] + 1;

//	Number of elements in banded global stiffness matrix
	NWK = DiagonalAddress[NEQ] - DiagonalAddress[0];

#ifdef _DEBUG_
	COutputter* Output = COutputter::Instance();
	//Output->PrintDiagonalAddress();
#endif

}

//	Assemble the banded gloabl stiffness matrix
void CDomain::AssembleStiffnessMatrix()
{
//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp.GetElement(0).SizeOfStiffnessMatrix();
		double* Matrix = new double[size];

//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
			ElementGrp.GetElement(Ele).assembly(Matrix, StiffnessMatrix);

		delete[] Matrix;
		Matrix = nullptr;
	}

#ifdef _DEBUG_
	COutputter* Output = COutputter::Instance();
	//Output->PrintStiffnessMatrix();
#endif

}

inline void sequeeze(sortnode*& list,int& len){
    QuickSort(list, 0, len);
    int lennew=1;
    for(unsigned int i=1;i<len;i++){if(list[i].key!=list[i-1].key)++lennew;}
    sortnode* newlist=new sortnode[lennew];
    newlist[0]=list[0];
    int iold=1,inew=0;
    while(iold<len){
        if(list[iold].key==newlist[inew].key){
            newlist[inew].value+=list[iold++].value;
        }
        else{
            newlist[++inew]=list[iold++];
        }
    }
    len=lennew;
    delete[] list;
    list=newlist;
    newlist=nullptr;
}


void CDomain::AssembleSparseSymmetricStiffnessMatrix(){
	//    Allocate for global force/displacement vector
    Force = new double[NEQ];
    clear(Force, NEQ);

    int* end_of_row=new int[NEQ];//at the first, it indices how many nnzs in each row, then we summerize them.
    clear(end_of_row,NEQ);
    //Symbol pre-indexing, which is useful for memory control;
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int size = ElementGrp.GetElement(0).ND;
#pragma parallel omp for
        for (unsigned int Ele = 0; Ele < ElementGrp.GetNUME(); Ele++){
			ElementGrp.GetElement(Ele).GenerateLocationMatrix();
            unsigned int* locationmatrix=ElementGrp.GetElement(Ele).LocationMatrix;
			QuickSort(locationmatrix, 0, size);
            for(int i=0;i<size;i++){
                if(locationmatrix[i])end_of_row[locationmatrix[i]-1]+=size-i;
            }
        }
	}
    
    sortnode** pnode=new sortnode*[NEQ];
#pragma parallel omp for
    for(int item=0;item<NEQ;item++){
        pnode[item]=new sortnode[end_of_row[item]];
    }
    
	int* index_of_row=new int[NEQ];
	for(int i=0;i<NEQ;i++){index_of_row[i]=0;}
	//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
		unsigned int size = ElementGrp.GetElement(0).ND;
		//	Loop over for all elements in group EleGrp
#pragma parallel omp for
        
		for (unsigned int Ele = 0; Ele < ElementGrp.GetNUME(); Ele++){
            double* Matrix=new double [300];
			ElementGrp.GetElement(Ele).GenerateLocationMatrix();
			unsigned int* locationmatrix=ElementGrp.GetElement(Ele).LocationMatrix;
            
			ElementGrp.GetElement(Ele).ElementStiffness(Matrix);
            
            for(int dy_index=0;dy_index<((size*(size+1))>>1);dy_index++){
                int i=0;
                while(((i*(i+1))>>1)<=dy_index){++i;}--i;
                int j=((i*(i+3))>>1)-dy_index;

                if(locationmatrix[i]>locationmatrix[j]){
                    if(!locationmatrix[j]){continue;}
                    int item=locationmatrix[j]-1;

                    pnode[item][index_of_row[item]].key=locationmatrix[i]-1;
                    pnode[item][index_of_row[item]++].value=Matrix[dy_index];
                    
                }
                else{
                    if(!locationmatrix[i]){continue;}
                    int item=locationmatrix[i]-1;
                    
                    pnode[item][index_of_row[item]].key=locationmatrix[j]-1;
                    pnode[item][index_of_row[item]++].value=Matrix[dy_index];

                }
                    //index 为假设用满对称矩阵存储时的位置。按行首尾相接。
            }
            
            delete[] Matrix;
		}
	}

#pragma parallel for
    for(int i=0;i<NEQ;i++){
        sequeeze(pnode[i], end_of_row[i]);
        //cout <<group_num[i]<< endl;
    }
 
	int* ia=new int[NEQ+1]();
    ia[0]=0;
    for(int i=0;i<NEQ;i++){ia[i+1]=ia[i]+end_of_row[i];}
    
    double* Matrix=new double[ia[NEQ]];
    int * ja=new int[ia[NEQ]];

	int inew=0;//用来遍历紧缩前和紧缩后的数组
    for(int i=0;i<NEQ;i++){
        for(int j=0;j<end_of_row[i];j++){
            Matrix[inew]=pnode[i][j].value;
            ja[inew++]=pnode[i][j].key;
        }
        delete[] pnode[i];
    }
    
    delete[] pnode; pnode=nullptr;
    delete[] end_of_row; end_of_row=nullptr;
	delete[] index_of_row; index_of_row=nullptr;

	solver=new Solver_Sparse(NEQ,ia,ja,Matrix);

}


//	Assemble the global nodal force vector for load case LoadCase
bool CDomain::AssembleForce(unsigned int LoadCase)
{
	if (LoadCase > NLCASE) 
		return false;

	CLoadCaseData* LoadData = &LoadCases[LoadCase - 1];

//	Loop over for all concentrated loads in load case LoadCase
	for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)
	{
		unsigned int dof = NodeList[LoadData->node[lnum] - 1].bcode[LoadData->dof[lnum] - 1];
        
        if(dof) // The DOF is activated
            Force[dof - 1] += LoadData->load[lnum];
	}
//    add the gravity
    for(unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
	    CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp.GetElement(0).ND;
		double* bodyforce = new double[size];

//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
			ElementGrp.GetElement(Ele).assemblygravity(bodyforce, Force, Gravity);

		delete[] bodyforce;
		bodyforce = nullptr;
    }
	return true;
}

//  Halfband optimization using GPS method
void CDomain::GPS()
{
	CBandwidth gps;
	gps.InitialBand();
// calculate the group
	gps.CalculateGraph();
// find the diameter points of the graph
	gps.Diameter();
// put nodes in layers
	gps.Layer();
// re numbering
	Order = new unsigned int[NEQ];
	gps.Renumber(Order);
	RefreshEquationNumber();
}

//	Allocate storage for matrices Force, ColumnHeights, DiagonalAddress and StiffnessMatrix
//	and calculate the column heights and address of diagonal elements
void CDomain::AllocateMatrices()
{
//	Allocate for global force/displacement vector
	Force = new double[NEQ];
    clear(Force, NEQ);

//  Create the banded stiffness matrix
    StiffnessMatrix = new CSkylineMatrix<double>(NEQ);

//	Calculate column heights
	CalculateColumnHeights();

//	Calculate address of diagonal elements in banded matrix
	CalculateDiagnoalAddress();

//	Allocate for banded global stiffness matrix
    StiffnessMatrix->Allocate();

}

