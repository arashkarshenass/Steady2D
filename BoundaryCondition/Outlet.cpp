/*=================================================
                Outlet Boundary Class source
=================================================*/

#include "Outlet.h"
#include "../Inputs.h"
#include "../Solver/OneCell.h"

#include <cmath>
#include <iostream>

using namespace std;

void Outlet::Initializer(double* wfc,double* wcn)
{
	wFaceCurr=wfc;
	wGhost=wcn;
}


void Outlet::CalculateGhost()
{
	double roD=wFaceCurr[0];
	double uD=abs(wFaceCurr[1]/roD);
	double pD=wFaceCurr[4];
	double aD=sqrt(gamma*pD/roD);

	//predictor
	double pG=pOut;
	double roG=roD+(pG-pD)/aD/aD;
	double uG=abs(uD-(pG-pD)/roD/aD);
	double aG=sqrt(gamma*pG/roG);

	//corrector
	double roAvg=(roG+roD)/2;
	double aAvg=(aG+aD)/2;
	roG=roD+(pG-pD)/aAvg/aAvg;
	uG=abs(uD-(pG-pD)/roAvg/aAvg);

	wGhost[0]=roG;
	wGhost[1]=roG*uG;
	wGhost[2]=0;
	wGhost[3]=pG/(gamma-1)+0.5*roG*uG*uG;
	wGhost[4]=pG;
}


