/*=================================================
                AUSM Class source
=================================================*/
#include "AUSM.h"
#include "../Inputs.h"

#include<iostream>
#include<cmath>
#include<algorithm>
using namespace std;



void AUSM::Initializer(double* wfc,double* wfn,double* ff,double**DL)
{
	wFaceCurr=wfc;
	wFaceNeigh=wfn;
	faceFlux=ff;
	dL=DL;
	fluxP=new double[4];
	fluxM=new double[4];
}

void AUSM::CalculateFaceFlux(int cNum,int fNum)
{
	//--------------current side
	ro=wFaceCurr[0];
	u=wFaceCurr[1]/ro;
	v=wFaceCurr[2]/ro;
	p=wFaceCurr[4];
	a=sqrt(gamma*p/ro);
	m=u/a;

	if(m<-1)
	{
		mP=0;
		pP=0;
	}
	else if(m>1)
	{
		mP=m;
		pP=p;
	}
	else
	{
		mP=((m+1)/2)*((m+1)/2);
		pP=p/2*(1+m);
	}

	fluxP[0]=mP*a*ro;
	fluxP[1]=mP*a*ro*u+pP;
	fluxP[2]=mP*a*ro*v;
	fluxP[3]=mP*a*(wFaceCurr[3]+p);

	//--------------neighbor side
	ro=wFaceNeigh[0];
	u=wFaceNeigh[1]/ro;
	v=wFaceNeigh[2]/ro;
	p=wFaceNeigh[4];
	a=sqrt(gamma*p/ro);
	m=u/a;

	if(m<-1)
	{
		mM=m;
		pM=p;
	}
	else if(m>1)
	{
		mM=0;
		pM=0;
	}
	else
	{
		mM=-((m-1)/2)*((m-1)/2);
		pM=p/2*(1-m);
	}

	fluxM[0]=mM*a*ro;
	fluxM[1]=mM*a*ro*u+pM;
	fluxM[2]=mM*a*ro*v;
	fluxM[3]=mM*a*(wFaceNeigh[3]+p);

	//--------------summing F+ and F--------------------
	faceFlux[0]=(fluxP[0]+fluxM[0])*dL[cNum][fNum];
	faceFlux[1]=(fluxP[1]+fluxM[1])*dL[cNum][fNum];
	faceFlux[2]=(fluxP[2]+fluxM[2])*dL[cNum][fNum];
	faceFlux[3]=(fluxP[3]+fluxM[3])*dL[cNum][fNum];

}
