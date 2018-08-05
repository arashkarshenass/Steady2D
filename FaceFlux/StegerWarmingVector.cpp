/*=================================================
                Steger Class source
=================================================*/
#include "../FaceFlux/StegerWarmingVector.h"
#include "../Inputs.h"

#include<iostream>
#include<cmath>
#include<algorithm>
using namespace std;

void StegerWarmingVector::Initializer(double* wfc,double* wfn,double* ff,double**DL)
{
	wFaceCurr=wfc;
	wFaceNeigh=wfn;
	faceFlux=ff;
	dL=DL;
	g1=gamma-1;
}

void StegerWarmingVector::CalculateFaceFlux(int cellNum,int faceNum)
{
	double rho,un,ut,a,s2;
	double lam[4]={0};
	double fluxPlus[4]={0};
	double fluxMinus[4]={0};

	//F+ (flux from current cell in n-t coord)
	rho=wFaceCurr[0];					//density
	un=wFaceCurr[1]/rho;				//u is Un (normal velocity)
	ut=wFaceCurr[2]/rho;				//v is Ut (tangential velocity)
	s2=un*un+ut*ut;						//square of speed
	a=sqrt(gamma*wFaceCurr[4]/rho);	//sound speed

	lam[0]=0.5*(un+abs(un));
	lam[1]=lam[0];
	lam[2]=0.5*((un+a)+abs(un+a));
	lam[3]=0.5*((un-a)+abs(un-a));

	fluxPlus[0]=rho/2/gamma*(2*g1*lam[1]+lam[2]+lam[3]);
	fluxPlus[1]=rho/2/gamma*(2*un*g1*lam[1]+(un+a)*lam[2]+(un-a)*lam[3]);
	fluxPlus[2]=rho/2/gamma*(2*ut*g1*lam[0]+ut*(lam[2]+lam[3]));
	fluxPlus[3]=rho/2/gamma*(2*ut*ut*g1*lam[0] + (un*un-ut*ut)*g1*lam[1] + (s2/2+a*un+a*a/g1)*lam[2] + (s2/2-a*un+a*a/g1)*lam[3]);

	//F- (flux from neighbor cell in n-t coord)
	rho=wFaceNeigh[0];
	un=wFaceNeigh[1]/rho;
	ut=wFaceNeigh[2]/rho;
	s2=un*un+ut*ut;
	a=sqrt(gamma*wFaceNeigh[4]/rho);

	lam[0]=0.5*(un-abs(un));
	lam[1]=lam[0];
	lam[2]=0.5*((un+a)-abs(un+a));
	lam[3]=0.5*((un-a)-abs(un-a));

	fluxMinus[0]=rho/2/gamma*(2*g1*lam[1]+lam[2]+lam[3]);
	fluxMinus[1]=rho/2/gamma*(2*un*g1*lam[1]+(un+a)*lam[2]+(un-a)*lam[3]);
	fluxMinus[2]=rho/2/gamma*(2*ut*g1*lam[0]+ut*(lam[2]+lam[3]));
	fluxMinus[3]=rho/2/gamma*(2*ut*ut*g1*lam[0] + (un*un-ut*ut)*g1*lam[1] + (s2/2+a*un+a*a/g1)*lam[2] + (s2/2-a*un+a*a/g1)*lam[3]);

	//summing F+ and F-
	faceFlux[0]=(fluxPlus[0]+fluxMinus[0])*dL[cellNum][faceNum];
	faceFlux[1]=(fluxPlus[1]+fluxMinus[1])*dL[cellNum][faceNum];
	faceFlux[2]=(fluxPlus[2]+fluxMinus[2])*dL[cellNum][faceNum];
	faceFlux[3]=(fluxPlus[3]+fluxMinus[3])*dL[cellNum][faceNum];
}

