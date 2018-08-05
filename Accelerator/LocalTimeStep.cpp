/*=================================================
                     LocalTimeStep class source
=================================================*/

#include "LocalTimeStep.h"
#include "../InputOutput/Mesh.h"
#include "../Inputs.h"

#include <iostream>
#include <cmath>
#include <iomanip>


using namespace std;

LocalTimeStep::LocalTimeStep(double** wTemp,Mesh* mp, double* dtTemp,int n)
{
	w=wTemp;								//variables coming from solver
	minDis=mp->GetMinDis();					//get minimum distance of each cell  from mesh
	dt=dtTemp;								//local time step vector coming from solver
	celltot=n;								//total cell number coming from solver
}

void LocalTimeStep::ClaculateLTS()
{
	u=0; v=0; flowSpeed=0; soundSpeed=0; totalSpeed=0;
	for(int i=0;i<celltot;i++)
	{
		u=w[i][1]/w[i][0];
		v=w[i][2]/w[i][0];
		flowSpeed=sqrt(u*u+v*v);						//calculating absolute flow speed of each cell
		soundSpeed=sqrt(gamma*w[i][4]/w[i][0]);			//calculating sound speed of each cell
		totalSpeed=flowSpeed+soundSpeed;				//calculating absolute flow speed + sound speed of each cell
		dt[i]=cfl*minDis[i]/totalSpeed;					//calculating local time step of each cell
	}
}



