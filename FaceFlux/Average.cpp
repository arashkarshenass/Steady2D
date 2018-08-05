/*=================================================
                Average Class source
=================================================*/
#include "Average.h"
#include "../Inputs.h"

#include<iostream>
#include<cmath>
#include<algorithm>
using namespace std;

void Average::Initializer(double* wfc,double* wfn,double* ff,double**DL)
{
	wFaceCurr=wfc;
	wFaceNeigh=wfn;
	faceFlux=ff;
	dL=DL;
}

void Average::CalculateFaceFlux(int cNum,int fNum)
{
	double ro=(wFaceCurr[0]+wFaceNeigh[0])/2;
	double u=(wFaceCurr[1]/wFaceCurr[0]+wFaceNeigh[1]/wFaceNeigh[0])/2;
	double v=(wFaceCurr[2]/wFaceCurr[0]+wFaceNeigh[2]/wFaceNeigh[0])/2;
	double etot=(wFaceCurr[3]/wFaceCurr[0]+wFaceNeigh[3]/wFaceNeigh[0])/2;
	double p=(wFaceCurr[4]+wFaceNeigh[4])/2;

	faceFlux[0]=ro*u*dL[cNum][fNum];
	faceFlux[1]=(ro*u*u+p)*dL[cNum][fNum];
	faceFlux[3]=ro*u*v*dL[cNum][fNum];
	faceFlux[4]=(ro*etot+p)*u*dL[cNum][fNum];
}

