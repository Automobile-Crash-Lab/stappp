/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <string>
#include <iostream>

#include "Domain.h"
#include "Bar.h"
#include "4Q.h"
#include "3T.h"
#include "8H.h"
#include "Beam.h"
#include "Plate.h"
#include "Shell.h"
#include "Outputter.h"
#include "Clock.h"

using namespace std;

int main(int argc, char *argv[])
{
	if (argc != 2) //  Print help message
	{
	    cout << "Usage: stap++ InputFileName\n";
		exit(1);
	}

	string filename(argv[1]);
    size_t found = filename.find_last_of('.');

    // If the input file name is provided with an extension
    if (found != std::string::npos) {
        if (filename.substr(found) == ".dat")
            filename = filename.substr(0, found);
        else {
            // The input file name must has an extension of 'dat'
            cout << "*** Error *** Invalid file extension: "
                 << filename.substr(found+1) << endl;
            exit(1);
        }
    }

    string InFile = filename + ".dat";
	string OutFile = filename + ".out";

	CDomain* FEMData = CDomain::Instance();

    Clock timer;
    timer.Start();

//  Read data and define the problem domain
	if (!FEMData->ReadData(InFile, OutFile))
	{
		cerr << "*** Error *** Data input failed!" << endl;
		exit(1);
	}
    
    double time_input = timer.ElapsedTime();

//  Bandwidth optimization using GPS method
//  Bandwidth optimization will only be operated if MODEX == 2
	if(FEMData->GetMODEX() == 2)
		FEMData->GPS();

//  Allocate global vectors and matrices, such as the Force, ColumnHeights,
//  DiagonalAddress and StiffnessMatrix, and calculate the column heights
//  and address of diagonal elements
    double time_assemble = timer.ElapsedTime();
    double* time_solution;
    double* time_stress;
    time_solution = new double[FEMData->GetNLCASE()];
    time_stress = new double[FEMData->GetNLCASE()];
    COutputter* Output = COutputter::Instance();
    
    
//  Assemble the banded gloabl stiffness matrix
if(FEMData->GetMODEX() == 3){
	FEMData->AssembleSparseSymmetricStiffnessMatrix();
    time_assemble = timer.ElapsedTime();
    Solver_Sparse* Solver_SP=FEMData->GetSparseSolver();
    Solver_SP->LDLT();

    //  Loop over for all load cases
    for (unsigned int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)
    {
        // Assemble righ-hand-side vector (force vector)
        FEMData->AssembleForce(lcase + 1);
        
        Solver_SP->BackSubstitution(FEMData->GetForce());
        
        time_solution[lcase] = timer.ElapsedTime();
        
        //Output->PrintDisplacement(lcase);
        
        Output->OutputNodalDisplacement(lcase);
    
        
        // Calculate and output stresses of all elements
        //Output->OutputElementStress();
        
        time_stress[lcase] = timer.ElapsedTime();
    }
    
    timer.Stop();
    
    //释放内存空间
    Solver_SP->release();
    
}
else{
    FEMData->AllocateMatrices();
    FEMData->AssembleStiffnessMatrix();

//  Solve the linear equilibrium equations for displacements
    
	CLDLTSolver* Solver = new CLDLTSolver(FEMData->GetStiffnessMatrix());
    
//  Perform L*D*L(T) factorization of stiffness matrix
    Solver->LDLT();

#ifdef _DEBUG_
    //Output->PrintStiffnessMatrix();
#endif

//  Loop over for all load cases
    for (unsigned int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)
    {
        // Assemble righ-hand-side vector (force vector)
        FEMData->AssembleForce(lcase + 1);

        //  Reduce right-hand-side force vector and back substitute
        Solver->BackSubstitution(FEMData->GetForce());

#ifdef _DEBUG_
        Output->PrintDisplacement(lcase);
#endif
            
        Output->OutputNodalDisplacement(lcase);

        time_solution[lcase] = timer.ElapsedTime();

        // Calculate and output stresses of all elements
	    Output->OutputElementStress();

        time_stress[lcase] = timer.ElapsedTime();
    }
    
    timer.Stop();
}
    
    *Output << "\n S O L U T I O N   T I M E   L O G   I N   S E C \n\n"
            << "     TIME FOR INPUT PHASE = " << time_input << endl
            << "     TIME FOR CALCULATION OF STIFFNESS MATRIX = " << time_assemble - time_input << endl;
    double time_solution_total = 0;
    for(unsigned int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)
    {
        double time_lcase;
        if (lcase == 0)
        {
            time_lcase = time_solution[lcase] - time_assemble;
            *Output << "     TIME FOR FACTORIZATION AND LOAD CASE SOLUTIONS" << endl
            << "     LOADCASE = " << lcase + 1 << "       TIME = " << time_lcase << endl;
        }
        else
        {
            time_lcase = time_solution[lcase] - time_solution[lcase-1];
            *Output << "     LOADCASE = " << lcase + 1 << "       TIME = " << time_lcase << endl;
        }
        time_solution_total += time_lcase;
    }
            *Output << "     TOTAL TIME FOR FACTORIZATION AND LOAD CASE SOLUTIONS = " << time_solution_total << endl << endl
            << "     T O T A L   S O L U T I O N   T I M E = " << time_stress[FEMData->GetNLCASE() - 1] << endl;

	return 0;
}
