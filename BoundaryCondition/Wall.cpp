/*=================================================
                Wall Boundary Class source
=================================================*/
#include "Wall.h"

void Wall::Initializer(double* wfc,double* wfn)
{
	wFaceCurr=wfc;
	wFaceNeigh=wfn;
}

void Wall::CalculateGhost()
{
	wFaceNeigh[0]=wFaceCurr[0];
	wFaceNeigh[1]=-wFaceCurr[1];
	wFaceNeigh[2]=wFaceCurr[2];
	wFaceNeigh[3]=wFaceCurr[3];
	wFaceNeigh[4]=wFaceCurr[4];
}

