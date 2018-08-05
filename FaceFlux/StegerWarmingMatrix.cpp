#include "StegerWarmingMatrix.h"
#include "../Inputs.h"

#include<iostream>
#include<cmath>
#include<algorithm>
using namespace std;

void StegerWarmingMatrix::Initializer(double*wfc,double*wfn,double*ff,double**dlT){
	wFaceCurr=wfc;
	wFaceNeigh=wfn;
	faceFlux=ff;
	dL=dlT;
	lam=new double[4];
	jacob=new double*[4];
	for(int i=0;i<4;i++)
		jacob[i]=new double[4];
	g1=gamma-1;

}

void StegerWarmingMatrix::CalculateFaceFlux(int cn,int fn){

	double fCurr[4]={0};
	double fNeigh[4]={0};

	//F from current cell
	un=wFaceCurr[1]/wFaceCurr[0];				ut=wFaceCurr[2]/wFaceCurr[0];
	a=sqrt(gamma*wFaceCurr[3]/wFaceCurr[0]);	s2=un*un+ut*ut;

	lam[0]=0.5*(un+abs(un));
	lam[1]=lam[0];
	lam[2]=0.5*((un+a)+abs(un+a));
	lam[3]=0.5*((un-a)+abs(un-a));

	CalculateacobianMatrix();
	fCurr[0]=jacob[0][0]*wFaceCurr[0] + jacob[0][1]*wFaceCurr[1] + jacob[0][2]*wFaceCurr[2] + jacob[0][3]*wFaceCurr[3];
	fCurr[1]=jacob[1][0]*wFaceCurr[0] + jacob[1][1]*wFaceCurr[1] + jacob[1][2]*wFaceCurr[2] + jacob[1][3]*wFaceCurr[3];
	fCurr[2]=jacob[2][0]*wFaceCurr[0] + jacob[2][1]*wFaceCurr[1] + jacob[2][2]*wFaceCurr[2] + jacob[2][3]*wFaceCurr[3];
	fCurr[3]=jacob[3][0]*wFaceCurr[0] + jacob[3][1]*wFaceCurr[1] + jacob[3][2]*wFaceCurr[2] + jacob[3][3]*wFaceCurr[3];

	//F from neighbor cell
	un=wFaceNeigh[1]/wFaceNeigh[0];				ut=wFaceNeigh[2]/wFaceNeigh[0];
	a=sqrt(gamma*wFaceNeigh[3]/wFaceNeigh[0]);	s2=un*un+ut*ut;

	lam[0]=0.5*(un-abs(un));
	lam[1]=lam[0];
	lam[2]=0.5*((un+a)-abs(un+a));
	lam[3]=0.5*((un-a)-abs(un-a));

	CalculateacobianMatrix();
	fNeigh[0]=jacob[0][0]*wFaceNeigh[0] + jacob[0][1]*wFaceNeigh[1] + jacob[0][2]*wFaceNeigh[2] + jacob[0][3]*wFaceNeigh[3];
	fNeigh[1]=jacob[1][0]*wFaceNeigh[0] + jacob[1][1]*wFaceNeigh[1] + jacob[1][2]*wFaceNeigh[2] + jacob[1][3]*wFaceNeigh[3];
	fNeigh[2]=jacob[2][0]*wFaceNeigh[0] + jacob[2][1]*wFaceNeigh[1] + jacob[2][2]*wFaceNeigh[2] + jacob[2][3]*wFaceNeigh[3];
	fNeigh[3]=jacob[3][0]*wFaceNeigh[0] + jacob[3][1]*wFaceNeigh[1] + jacob[3][2]*wFaceNeigh[2] + jacob[3][3]*wFaceNeigh[3];

	//face flux
	faceFlux[0]=(fCurr[0]+fNeigh[0])*dL[cn][fn];
	faceFlux[1]=(fCurr[1]+fNeigh[1])*dL[cn][fn];
	faceFlux[2]=(fCurr[2]+fNeigh[2])*dL[cn][fn];
	faceFlux[3]=(fCurr[3]+fNeigh[3])*dL[cn][fn];
}

void StegerWarmingMatrix::CalculateacobianMatrix(){
	jacob[0][0]=lam[0]* (1-g1*s2/2/a/a) + lam[2]* ((g1*s2-2*un*a)/4/a/a) + lam[3]* ((g1*s2+2*un*a)/4/a/a);
	jacob[0][1]=lam[0]* (g1*un/a/a) + lam[2]* ((a-g1*un)/2/a/a) + lam[3]* (-(a+g1*un)/2/a/a);
	jacob[0][2]=g1*ut/2/a/a* (2*lam[0] - lam[2] - lam[3]);
	jacob[0][3]=-g1/2/a/a* (2*lam[0] - lam[2] - lam[3]);
	jacob[1][0]=lam[0]* ((2*a*a-g1*s2)*un/2/a/a) + lam[2]* ((un+a)*(g1*s2-2*a*un)/4/a/a) + lam[3]* ((un-a)*(g1*s2+2*a*un)/4/a/a);
	jacob[1][1]=lam[0]* (g1*un*un/a/a) + lam[2]* ((a+un)*(a-g1*un)/2/a/a) + lam[3]* ((a-un)*(a+g1*un)/2/a/a);
	jacob[1][2]=lam[0]* (g1*un*ut/a/a) + lam[2]* (-g1*(un+a)*ut/2/a/a) + lam[3]* (-g1*(un-a)*ut/2/a/a);
	jacob[1][3]=lam[0]* (-g1*un/a/a)+ lam[2]* (g1*(un+a)/2/a/a) + lam[3]* (g1*(un-a)/2/a/a);
	jacob[2][0]=lam[0]* (-g1*ut*s2/2/a/a) + lam[2]* ((g1*s2-2*un*a)*ut/4/a/a) + lam[3]*  ((g1*s2+2*un*a)*ut/4/a/a);
	jacob[2][1]=lam[0]* (g1*un*ut/a/a) + lam[2]* ((a-g1*un)*ut/2/a/a) + lam[3]* (-(a+g1*un)*ut/2/a/a) ;
	jacob[2][2]=g1*ut*ut/2/a/a* (2*lam[0] - lam[2] - lam[3]) + lam[0];
	jacob[2][3]=-g1*ut/2/a/a* (2*lam[0] - lam[2] - lam[3]);
	jacob[3][0]=lam[0]* ((-g1*s2*s2+2*a*a*(un*un-ut*ut))/4/a/a) + lam[2]* ((g1*s2-2*a*un)*(2*a*a/g1+s2+2*a*un)/8/a/a) + lam[3]* ((g1*s2+2*a*un)*(2*a*a/g1+s2-2*a*un)/8/a/a);
	jacob[3][1]=lam[0]* (g1*un*s2/2/a/a) + lam[2]* ((2*a*a/g1+s2+2*a*un)*(a-g1*un)/4/a/a) + lam[3]* (-(2*a*a/g1+s2-2*a*un)*(a+g1*un)/4/a/a);
	jacob[3][2]=lam[0]* ((g1*s2+2*a*a)*ut/2/a/a) + lam[2]* (-ut*g1*(2*a*a/g1+s2+2*un*a)/4/a/a) + lam[3]* (-ut*g1*(2*a*a/g1+s2-2*un*a)/4/a/a);
	jacob[3][3]=lam[0]* (-g1*s2/2/a/a) + lam[2]* (g1*(2*a*a/g1+s2+2*un*a)/4/a/a) + lam[3]* (g1*(2*a*a/g1+s2-2*un*a)/4/a/a);
}
