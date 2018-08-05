/*=================================================
                Inlet Boundary Class source
=================================================*/

#include "Inlet.h"
#include "../Inputs.h"
#include "../Solver/OneCell.h"

#include <cmath>
#include <iostream>
using namespace std;

void Inlet::Initializer(double* wfc,double* wcn)
{
	wFaceCurr=wfc;
	wCellNeigh=wcn;
}


void Inlet::CalculateGhost()
{
	double erRef=0.0001;
	double roD=wFaceCurr[0];
	double uD=abs(wFaceCurr[1]/roD);
	double pD=wFaceCurr[4];
	double aD=sqrt(gamma*pD/roD);
	double y=0;
	double yp=0;

	//predictor
	double A=pD/roD/aD-uD;
	double B=pStagIn/roD/aD;
	double C=(gamma-1)/2/gamma/R/tStagIn;
	double phi=pD/pStagIn;

	double errorPhi=1;
	double phiOld=0;
	while(errorPhi>erRef){
		phiOld=phi;
		y=pow(phi,(gamma-1)/gamma)+B*B*C*phi*phi-2*A*B*C*phi+A*A*C-1;
		yp=(gamma-1)/gamma*pow(phi,-1/gamma)+2*B*B*C*phi-2*A*B*C;
		phi=phi-y/yp;
		errorPhi=abs((phi-phiOld)/phi);
	}

	double pG=pStagIn*phi;
	double roG=pStagIn/R/tStagIn*pow(phi,1/gamma);
	double uG=abs(-A+B*phi);
	double aG=sqrt(gamma*pG/roG);

	//corrector
	double roAvg=(roG+roD)/2;
	double aAvg=(aG+aD)/2;
	A=pD/roAvg/aAvg-uD;
	B=pStagIn/roAvg/aAvg;

	errorPhi=1;
	phiOld=0;
	while(errorPhi>erRef){
		phiOld=phi;
		y=pow(phi,(gamma-1)/gamma)+B*B*C*phi*phi-2*A*B*C*phi+A*A*C-1;
		yp=(gamma-1)/gamma*pow(phi,-1/gamma)+2*B*B*C*phi-2*A*B*C;
		phi=phi-y/yp;
		errorPhi=abs((phi-phiOld)/phi);
	}

	pG=pStagIn*phi;
	roG=pStagIn/R/tStagIn*pow(phi,1/gamma);
	uG=abs(-A+B*phi);

/*	double errorU=1;
	double uOld=0;
	while(errorU>erRef)
	{
		uOld=uG;
		roAvg=(roG+roD)/2;
		aAvg=(aG+aD)/2;
		A=pD/roAvg/aAvg-uD;
		B=pStagIn/roAvg/aAvg;

		errorPhi=1;
		phiOld=0;
		while(errorPhi>erRef){
			phiOld=phi;
			y=pow(phi,(gamma-1)/gamma)+B*B*C*phi*phi-2*A*B*C*phi+A*A*C-1;
			yp=(gamma-1)/gamma*pow(phi,-1/gamma)+2*B*B*C*phi-2*A*B*C;
			phi=phi-y/yp;
			errorPhi=abs((phi-phiOld)/phi);
		}

		pG=pStagIn*phi;
		roG=pStagIn/R/tStagIn*pow(phi,1/gamma);
		uG=abs(-A+B*phi);
		errorU=abs((uG-uOld)/uG);
	}*/

	wCellNeigh[0]=roG;
	wCellNeigh[1]=roG*uG;
	wCellNeigh[2]=0;
	wCellNeigh[3]=pG/(gamma-1)+0.5*roG*uG*uG;
	wCellNeigh[4]=pG;
}


