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
#include "Outputter.h"
#include "SkylineMatrix.h"

#include <iostream>
#include <iomanip>
#include <ctime>

using namespace std;

//	Output current time and date
void COutputter::PrintTime(const struct tm* ptm, COutputter &output)
{
	const char* weekday[] = {"Sunday", "Monday", "Tuesday", "Wednesday",
							 "Thursday", "Friday", "Saturday"};
	const char* month[] = {"January", "February", "March", "April", "May", "June",
						   "July", "August", "September", "October", "November", "December"};

	output << "        (";
	output << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec << " on ";
	output << month[ptm->tm_mon] << " " << ptm->tm_mday << ", " << ptm->tm_year + 1900 << ", "
		   << weekday[ptm->tm_wday] << ")" << endl
		   << endl;
}

COutputter* COutputter::_instance = nullptr;

//	Constructor
COutputter::COutputter(string FileName)
{
	OutputFile.open(FileName);

	if (!OutputFile)
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}
}

//	Return the single instance of the class
COutputter* COutputter::Instance(string FileName)
{
	if (!_instance)
		_instance = new COutputter(FileName);
	return _instance;
}

//	Print program logo
void COutputter::OutputHeading()
{
	CDomain* FEMData = CDomain::Instance();

	*this << "TITLE : " << FEMData->GetTitle() << endl;

	time_t rawtime;
	struct tm* timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	PrintTime(timeinfo, *this);
}

//	Print nodal data
void COutputter::OutputNodeInfo()
{
	CDomain* FEMData = CDomain::Instance();

	CNode* NodeList = FEMData->GetNodeList();

	*this << "C O N T R O L   I N F O R M A T I O N" << endl
		  << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	unsigned int NUMNP = FEMData->GetNUMNP();
	unsigned int NUMEG = FEMData->GetNUMEG();
	unsigned int NLCASE = FEMData->GetNLCASE();
	unsigned int MODEX = FEMData->GetMODEX();

	*this << "      NUMBER OF NODAL POINTS . . . . . . . . . . (NUMNP)  =" << setw(6) << NUMNP << endl;
	*this << "      NUMBER OF ELEMENT GROUPS . . . . . . . . . (NUMEG)  =" << setw(6) << NUMEG << endl;
	*this << "      NUMBER OF LOAD CASES . . . . . . . . . . . (NLCASE) =" << setw(6) << NLCASE << endl;
	*this << "      SOLUTION MODE  . . . . . . . . . . . . . . (MODEX)  =" << setw(6) << MODEX << endl;
	*this << "         EQ.0, DATA CHECK" << endl
		  << "         EQ.1, EXECUTION" << endl
		  << endl;

	*this << " N O D A L   P O I N T   D A T A" << endl << endl;
	*this << "    NODE       BOUNDARY                         NODAL POINT" << endl
		  << "   NUMBER  CONDITION  CODES                     COORDINATES" << endl;

	for (unsigned int np = 0; np < NUMNP; np++)
		NodeList[np].Write(*this, np);

	*this << endl;
}

//	Output equation numbers
void COutputter::OutputEquationNumber()
{
	CDomain* FEMData = CDomain::Instance();
	unsigned int NUMNP = FEMData->GetNUMNP();

	CNode* NodeList = FEMData->GetNodeList();

	*this << " EQUATION NUMBERS" << endl
		  << endl;
	*this << "   NODE NUMBER   DEGREES OF FREEDOM" << endl;
	*this << "        N           X    Y    Z" << endl;

	for (unsigned int np = 0; np < NUMNP; np++) // Loop over for all node
		NodeList[np].WriteEquationNo(*this, np);

	*this << endl;
}

//	Output element data
void COutputter::OutputElementInfo()
{
//	Print element group control line

	CDomain* FEMData = CDomain::Instance();

	unsigned int NUMEG = FEMData->GetNUMEG();

	*this << " E L E M E N T   G R O U P   D A T A" << endl
		  << endl
		  << endl;

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		*this << " E L E M E N T   D E F I N I T I O N" << endl
			  << endl;

		ElementTypes ElementType = FEMData->GetEleGrpList()[EleGrp].GetElementType();
		unsigned int NUME = FEMData->GetEleGrpList()[EleGrp].GetNUME();

		*this << " ELEMENT TYPE  . . . . . . . . . . . . .( NPAR(1) ) . . =" << setw(5)
			  << ElementType << endl;
		*this << "     EQ.1, BAR ELEMENTS" << endl
			  << "     EQ.2, 4Q ELEMENTS" << endl
			  << "     EQ.3, 3T ELEMENTS" << endl
			  << "     EQ.4, 8H ELEMENTS" << endl
			  << "     EQ.5, BEAM ELEMENTS" << endl
			  << "     EQ.6, Plate ELEMENTS" << endl
			  << "     EQ.7, Shell ELEMENTS" << endl
			  << "     EQ.8, NOT AVILABLE CURRENTLY" << endl
			  << endl;

		*this << " NUMBER OF ELEMENTS. . . . . . . . . . .( NPAR(2) ) . . =" << setw(5) << NUME
			  << endl
			  << endl;

		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				PrintBarElementData(EleGrp);
				break;
			case ElementTypes::T3: // Bar element
				Print3TElementData(EleGrp);
				break;
			case ElementTypes::Q4: // 4Q element
				Print4QElementData(EleGrp);
				break;
			case ElementTypes::Plate: // 4Q element
				PrintPlateElementData(EleGrp);
				break;
            case ElementTypes::Shell: // 4Q element
                PrintShellElementData(EleGrp);
				break;
			case ElementTypes::H8:
				Print8HElementData(EleGrp);
				break;
			case ElementTypes::Beam: //Beam element
				PrintBeamElementData(EleGrp);
				break;
            default:
                break;
		}
	}
}
//	Output bar element data
void COutputter::PrintBarElementData(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::Instance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S     CROSS-SECTIONAL	  DENSITY" << endl
		  << " NUMBER     MODULUS          AREA" << endl
		  << "               E              A			  rho" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
		ElementGroup.GetMaterial(mset).Write(*this, mset);

	*this << endl
		  << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
	*this << " ELEMENT     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      I        J       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
		ElementGroup.GetElement(Ele).Write(*this, Ele);

	*this << endl;
}
void COutputter::Print3TElementData(unsigned int EleGrp){
	CDomain* FEMData = CDomain::Instance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S        POSSION     THICKNESS" << endl
		  << " NUMBER     MODULUS         RAITO" << endl
		  << "               E              nu           t" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
		ElementGroup.GetMaterial(mset).Write(*this, mset);

	*this << endl
		  << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
	*this << " ELEMENT     NODE     NODE     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      I        J        K        L       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
		ElementGroup.GetElement(Ele).Write(*this, Ele);

	*this << endl;
}

void COutputter::Print4QElementData(unsigned int EleGrp){
	CDomain* FEMData = CDomain::Instance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;
	
	*this << "  SET       YOUNG'S        POSSION     THICKNESS" << endl
		  << " NUMBER     MODULUS         RAITO" << endl
		  << "               E              nu           t" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
		ElementGroup.GetMaterial(mset).Write(*this, mset);

	*this << endl
		  << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
	*this << " ELEMENT     NODE     NODE     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      I        J        K        L       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
		ElementGroup.GetElement(Ele).Write(*this, Ele);

	*this << endl;
}
void COutputter::PrintPlateElementData(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::Instance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S        POSSION     THICKNESS" << endl
		  << " NUMBER     MODULUS         RAITO" << endl
		  << "               E              nu           t" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
		ElementGroup.GetMaterial(mset).Write(*this, mset);

	*this << endl
		  << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
	*this << " ELEMENT     NODE     NODE     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      I        J        K        L       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
		ElementGroup.GetElement(Ele).Write(*this, Ele);

	*this << endl;
}


//	Output beam element data
void COutputter::PrintBeamElementData(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::Instance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S      POISSON     CROSS-SECTIONAL     inertia         inertia          inertia" << endl
		 << " NUMBER     MODULUS        RATIO        AREA            moment_y        moment_z         moment_x" << endl
		 << "               E               v          A               Iy               Iz             Ix     " << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
		ElementGroup.GetMaterial(mset).Write(*this, mset);

	*this << endl
		  << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
	*this << " ELEMENT     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      I        J       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
		ElementGroup.GetElement(Ele).Write(*this, Ele);

	*this << endl;
}
void COutputter::PrintShellElementData(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::Instance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S        POSSION     THICKNESS" << endl
		  << " NUMBER     MODULUS         RAITO" << endl
		  << "               E              nu           t" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
		ElementGroup.GetMaterial(mset).Write(*this, mset);

	*this << endl
		  << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
	*this << " ELEMENT     NODE     NODE     NODE     NODE     NODE     NODE     NODE     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      1        2        3        4        5        6        7        8        9       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
		ElementGroup.GetElement(Ele).Write(*this, Ele);

	*this << endl;
}


//	Print load data
void COutputter::OutputLoadInfo()
{
	CDomain* FEMData = CDomain::Instance();

	for (unsigned int lcase = 1; lcase <= FEMData->GetNLCASE(); lcase++)
	{
		CLoadCaseData* LoadData = &FEMData->GetLoadCases()[lcase - 1];

		*this << setiosflags(ios::scientific);
		*this << " L O A D   C A S E   D A T A" << endl
			  << endl;

		*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
		*this << "     NUMBER OF CONCENTRATED LOADS . =" << setw(6) << LoadData->nloads << endl
			  << endl;
		*this << "    NODE       DIRECTION      LOAD" << endl
			  << "   NUMBER                   MAGNITUDE" << endl;

		LoadData->Write(*this, lcase);

		*this << endl;
	}
}

void COutputter::Print8HElementData(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::Instance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S         POSSION" << endl
		  << " NUMBER     MODULUS          RATIO" << endl
		  << "               E              NU" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
		ElementGroup.GetMaterial(mset).Write(*this, mset);

	*this << endl
		  << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
	*this << " ELEMENT     NODE     NODE     NODE     NODE     NODE     NODE     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      0        1        2        3        4        5        6        7       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
		ElementGroup.GetElement(Ele).Write(*this, Ele);

	*this << endl;
}


//	Print nodal displacement
void COutputter::OutputNodalDisplacement(unsigned int lcase)
{
	CDomain* FEMData = CDomain::Instance();
	CNode* NodeList = FEMData->GetNodeList();
	double* Displacement = FEMData->GetDisplacement();

	*this << " LOAD CASE" << setw(5) << lcase + 1 << endl
		  << endl
		  << endl;

	*this << setiosflags(ios::scientific);

	*this << " D I S P L A C E M E N T S" << endl
		  << endl;
	*this << "  NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT" << endl;

	for (unsigned int np = 0; np < FEMData->GetNUMNP(); np++)
		NodeList[np].WriteNodalDisplacement(*this, np, Displacement);

	*this << endl;
}

//	Calculate stresses
void COutputter::OutputElementStress()
{
	CDomain* FEMData = CDomain::Instance();

	double* Displacement = FEMData->GetDisplacement();

	unsigned int NUMEG = FEMData->GetNUMEG();

	for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
	{
		*this << " S T R E S S  C A L C U L A T I O N S  F O R  E L E M E N T  G R O U P" << setw(5)
			  << EleGrpIndex + 1 << endl
			  << endl;

		CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
		unsigned int NUME = EleGrp.GetNUME();
		ElementTypes ElementType = EleGrp.GetElementType();

		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				{*this << "  ELEMENT             FORCE            STRESS" << endl
					<< "  NUMBER" << endl;

				double stress;

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp.GetElement(Ele);
					Element.ElementStress(&stress, Displacement);

					CBarMaterial& material = *dynamic_cast<CBarMaterial*>(Element.GetElementMaterial());
					*this << setw(5) << Ele + 1 << setw(22) << stress * material.Area << setw(18)
						<< stress << endl;
				}

				*this << endl;

				break;
				}
				
			case ElementTypes::Q4:{
				*this << "  ELEMENT             FORCE            STRESS" << endl
					<< "  NUMBER" << endl;

				double* stress4Q=new double[4];
                
                for (unsigned int Ele = 0; Ele < NUME; Ele++)
                {
                    CElement& Element = EleGrp.GetElement(Ele);
                    Element.ElementStress(stress4Q, Displacement);
                    *this << setw(5) << Ele + 1;
                    for(int i=0;i<4;i++){
                        *this << setw(5) << (Element.GetNodes())[i]->NodeNumber << setw(14)
                        << stress4Q[i];
                    }
                }
                *this << endl;
                break;
            }

            case ElementTypes::T3:{
				*this << "  ELEMENT             FORCE            STRESS" << endl
					<< "  NUMBER" << endl;

				double* stress3T=new double[12];
                double zero=0.0;
                for (unsigned int Ele = 0; Ele < NUME; Ele++)
                {
                    C3T* Element = dynamic_cast<C3T*>(&(EleGrp.GetElement(Ele)));
                    Element->ElementStress(stress3T, Displacement);
                    *this << setw(5) << Ele + 1;
                    for(int ii=0;ii<4;ii++){
                    	*this << setw(5) << 1+ii<< setw(22) << Element->location[2*ii]
					<< setw(22) << Element->location[2*ii+1]<< setw(22)<< zero<<
					 setw(22) << stress3T[3*ii]
					<< setw(22) << stress3T[3*ii+1]<< setw(22) << stress3T[3*ii+2];
                    }
                }
                *this << endl;
                break;
            }
	        case ElementTypes::Beam: // beam element
				{*this << "  ELEMENT          STRESSx         Moment1         Moment2         Shear_force" << endl
					<< "  NUMBER" << endl;

				double* stressBeam=new double[4];

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp.GetElement(Ele);
					Element.ElementStress(stressBeam, Displacement);

					CBeamMaterial& material = *dynamic_cast<CBeamMaterial*>(Element.GetElementMaterial());
					*this << setw(5) << Ele + 1 << setw(22) << stressBeam[0] << setw(18) << stressBeam[1] << setw(18)
					<< stressBeam[2] << setw(18) << stressBeam[3] << setw(18)
					<< endl;
				}

				*this << endl;

				break;
                }
			case ElementTypes::Plate: // 4Q element
				{*this << "  ELEMENT          STRESSxx          STRESSyy          STRESSxy" << endl
					<< "  NUMBER" << endl;

				double stress[12];

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp.GetElement(Ele);
					Element.ElementStress(stress, Displacement);

					CPlateMaterial& material = *dynamic_cast<CPlateMaterial*>(Element.GetElementMaterial());
					*this << setw(5) << Ele + 1 << setw(22) << stress[0] << setw(18) << stress[1] << setw(18)<< stress[2] << setw(18) << endl;
					*this << setw(27) << stress[3] << setw(18) << stress[4] << setw(18) << stress[5] << endl;
					*this << setw(27) << stress[6] << setw(18) << stress[7] << setw(18) << stress[8] << endl;
					*this << setw(27) << stress[9] << setw(18) << stress[10] << setw(18) << stress[11] << endl;
				}

				*this << endl;

				break;}


			case ElementTypes::H8:
				*this << "  ELEMENT           STRESSxx           STRESSyy           STRESSzz           STRESSyz           STRESSzx           STRESSxy" << endl
					<< "  NUMBER" << endl;

				double stress4[6];

				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp.GetElement(Ele);
					Element.ElementStress(stress4, Displacement);

					C8HMaterial& material = *dynamic_cast<C8HMaterial*>(Element.GetElementMaterial());
					*this << setw(5) << Ele + 1 << setw(22)
						<< stress4[0] << setw(18)
						<< stress4[1] << setw(18)
						<< stress4[2] << setw(18)
						<< stress4[3] << setw(18)
						<< stress4[4] << setw(18)
						<< stress4[5] << setw(18)<<endl;
				}

				*this << endl;

				break;

			default: // Invalid element type
				cerr << "*** Error *** Elment type " << ElementType
					<< " has not been implemented.\n\n";
		}
	}
}

//	Print total system data
void COutputter::OutputTotalSystemData()
{
	CDomain* FEMData = CDomain::Instance();

	*this << "	TOTAL SYSTEM DATA" << endl
		  << endl;

	*this << "     NUMBER OF EQUATIONS . . . . . . . . . . . . . .(NEQ) = " << FEMData->GetNEQ()
		  << endl
		  << "     NUMBER OF MATRIX ELEMENTS . . . . . . . . . . .(NWK) = " << FEMData->GetNWK()
		  << endl
		  << "     MAXIMUM HALF BANDWIDTH  . . . . . . . . . . . .(MK ) = " << FEMData->GetMK()
		  << endl
		  << "     MEAN HALF BANDWIDTH . . . . . . . . . . . . . .(MM ) = "
		  << FEMData->GetNWK() / FEMData->GetNEQ() << endl
		  << endl
		  << endl;
}

#ifdef _DEBUG_

//	Print column heights for debuging
void COutputter::PrintColumnHeights()
{
	*this << "*** _Debug_ *** Column Heights" << endl;

	CDomain* FEMData = CDomain::Instance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* ColumnHeights = StiffnessMatrix->GetColumnHeights();

	for (unsigned int col = 0; col < NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << ColumnHeights[col];
	}

	*this << endl
		  << endl;
}

//	Print address of diagonal elements for debuging
void COutputter::PrintDiagonalAddress()
{
	*this << "*** _Debug_ *** Address of Diagonal Element" << endl;

	CDomain* FEMData = CDomain::Instance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	for (unsigned int col = 0; col <= NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << DiagonalAddress[col];
	}

	*this << endl
		  << endl;
}

//	Print banded and full stiffness matrix for debuging
void COutputter::PrintStiffnessMatrix()
{
	*this << "*** _Debug_ *** Banded stiffness matrix" << endl;

	CDomain* FEMData = CDomain::Instance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < DiagonalAddress[NEQ] - DiagonalAddress[0]; i++)
	{
		*this << setw(14) << (*StiffnessMatrix)(i);

		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}
	}

	*this << endl
		  << endl;

	*this << "*** _Debug_ *** Full stiffness matrix" << endl;

	for (int I = 1; I <= NEQ; I++)
	{
		for (int J = 1; J <= NEQ; J++)
		{
			int H = DiagonalAddress[J] - DiagonalAddress[J - 1];
			if (J - I - H >= 0)
			{
				*this << setw(14) << 0.0;
			}
			else
			{
				*this << setw(14) << (*StiffnessMatrix)(I, J);
			}
		}

		*this << endl;
	}

	*this << endl;
}

//	Print displacement vector for debuging
void COutputter::PrintDisplacement(unsigned int loadcase)
{
	*this << "*** _Debug_ *** Displacement vector" << endl;

	CDomain* FEMData = CDomain::Instance();

	unsigned int NEQ = FEMData->GetNEQ();
	double* Force = FEMData->GetForce();

	*this << "  Load case = " << loadcase << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < NEQ; i++)
	{
		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}

		*this << setw(14) << Force[i];
	}

	*this << endl
		  << endl;
}

#endif
