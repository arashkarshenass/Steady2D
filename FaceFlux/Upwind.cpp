/*=================================================
                Steger Class source
=================================================*/
#include "../FaceFlux/Upwind.h"
#include "../Inputs.h"

#include<iostream>
#include<cmath>
#include<algorithm>
using namespace std;

void Upwind::Initializer(double* wfc,double* wfn,double* ff,double** DL)
{
	wFaceCurr=wfc;
	wFaceNeigh=wfn;
	faceFlux=ff;
	dl=DL;
}

void Upwind::CalculateFaceFlux(int cNum,int fNum)
{
	double rhoCurr=wFaceCurr[0];
	double pCurr=wFaceCurr[4];
	double uCurr=wFaceCurr[1]/rhoCurr;
	double aCurr=sqrt(gamma*pCurr/rhoCurr);
	double mCurr=uCurr/aCurr;

	double rhoNeigh=wFaceNeigh[0];
	double pNeigh=wFaceNeigh[4];
	double uNeigh=wFaceNeigh[1]/rhoNeigh;
	double aNeigh=sqrt(gamma*pNeigh/rhoNeigh);
	double mNeigh=uNeigh/aNeigh;

	double mF=0.5*(mCurr+mNeigh);
	double rhoF=0;
	double etotF=0;
	double uF=0;
	double vF=0;
	double pF=0;
	double htotF=0;

	if (mF>0)			//comes from current cell
	{
		rhoF=rhoCurr;
		etotF=wFaceCurr[3]/rhoCurr;
		vF=wFaceCurr[2]/rhoCurr;
		if(mF>1)
		{
			pF=pCurr;
			uF=uCurr;
		}
		else
		{
			pF=(pCurr+pNeigh)/2;
			uF=(uCurr + uNeigh)/2;
		}
		htotF=etotF+(pF/rhoF);
	}
	else			//comes from Neigh cell
	{
		rhoF=rhoNeigh;
		etotF=wFaceNeigh[3]/rhoNeigh;
		vF=wFaceNeigh[2]/rhoNeigh;
		if(mF<-1)
		{
			pF=pNeigh;
			uF=uNeigh;
		}
		else
		{
			pF=(pCurr+pNeigh)/2;
			uF=(uCurr + uNeigh)/2;
		}
		htotF=etotF+(pF/rhoF);
	}

	faceFlux[0]=rhoF*uF*dl[cNum][fNum];
	faceFlux[1]=(rhoF*uF*uF+pF)*dl[cNum][fNum];
	faceFlux[2]=rhoF*uF*vF*dl[cNum][fNum];
	faceFlux[3]=rhoF*uF*htotF*dl[cNum][fNum];
}



