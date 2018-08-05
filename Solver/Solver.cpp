/*=================================================
                     solver class source
=================================================*/
#include "Solver.h"
#include "OneCell.h"
#include "Residue.h"
#include "../Inputs.h"
#include "../InputOutput/Mesh.h"
#include "../Accelerator/LocalTimeStep.h"
#include "../Utility/Utility.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

Solver::Solver(Mesh* mpTemp, Utility* upTemp)
{
	up=upTemp;												//pointer to utility object
	mp=mpTemp;												//pointer to mesh Object

	celltot=mp->GetCelltot();									//get total cell number from mesh object
	area=mp->GetArea();											//get area of each cell from mesh object

	w=new double*[celltot];										//allocating memory for conservative variables
	for (int i=0;i<celltot;i++)
		w[i]=new double[5];
	up->MatrixInitializer(w,celltot,5,0);							//set w=0

	wOld=new double*[celltot];									//allocating memory for old conservative variables
	for (int i=0;i<celltot;i++)
		wOld[i]=new double[5];
	up->MatrixInitializer(wOld,celltot,5,0);						//set wOld=0

	dw=new double*[celltot];									//allocating memory for conservative variables
	for (int i=0;i<celltot;i++)
		dw[i]=new double[4];
	up->MatrixInitializer(dw,celltot,4,0);							//set dw=0

	totalFlux=new double*[celltot];								//allocating memory for net fluxes of each cell (x-y coord)
	for (int i=0;i<celltot;i++)
		totalFlux[i]=new double[4];
	up->MatrixInitializer(totalFlux,celltot,4,0);					//set totalFlux=0

	dt=new double[celltot];										//allocating memory for local time steps
	up->VectorInitializer(dt,celltot,0);						//set dt=0

	res= new double[4];											//allocating memory for residue values of 4 equations
	up->VectorInitializer(res,4,0);								//set red=0

	//Initialization
	for (int i=0;i<celltot;i++)									//set initial value of coservative variables for all cells
	{	w[i][0]=ro0;											//inititial rho
		w[i][1]=ro0*u0;											//initial rho*u
		w[i][2]=ro0*v0;											//initial rho*v
		w[i][3]=p0/(gamma-1)+0.5*ro0*(u0*u0+v0*v0);				//initial rho*E
		w[i][4]=p0;												//initial P
	}

	if (prntsln==1){
	puts("\nInitialization is done as below for all cells:");
	puts(" ro      ro.u    	 ro.v       ro.e          P");
	puts("---------------------------------------------------");
	cout<<left<<setw(8)<<setprecision(2)<<fixed<<w[0][0]<<left<<setw(10)<<setprecision(2)<<fixed<<w[0][1] \
				<<left<<setw(13)<<setprecision(2)<<fixed<<w[0][2]<<left<<setw(12)<<setprecision(2)<<fixed<<w[0][3] \
				<<left<<setw(13)<<setprecision(2)<<fixed<<w[0][4]<<endl;
	}

	LocalTimeStep localTimeStep(w,mp,dt,celltot);			//create localTimeStep object with inputs to instructor
	tp=&localTimeStep;												//address of localTimeStep object

	OneCell oneCell(mp,w,totalFlux,up);							//create cellFlux object with inputs to its structor
	ocp=&oneCell;														//address of cellFlux object

	Residue residue(dw,res,up,celltot);							//create residue object and send addresses to it
	rp=&residue;															//address of residue object

	Iterator();
}


void Solver::Iterator()
{
	puts("\nIteration is started.");
	puts("\n-------------------Resuduals----------------");
	puts("#------------ro--------ro.u--------ro.v--------ro.e");
	double alpha[]={0.666,0.666,1};							//r-k3 coefficients
	itr=0;

	while (res[3]>resRef)									//psudo-time iteration
	{
		itr+=1;												//Iteration counter

		if(itr%nfdts==0 || itr==1)
			tp->ClaculateLTS();								//calculate local time step

		up->MatrixEqualizer(wOld,w,celltot,5);				//set wOld equal to w

		for (int rk=0;rk<3;rk++)								//r-k3 loop
		{
			up->MatrixInitializer(totalFlux,celltot,4,0);	//Initialize totalFlux to 0 at begining each rk iteration
			up->MatrixInitializer(dw,celltot,4,0);			//Initialize dw to 0 at begining each rk iteration

			for (int cellNum=0;cellNum<celltot;cellNum++){				//cell loop
				ocp->ClaculateCellFlux(cellNum);				//calculatig total fluxes of cell n
				for (int j=0;j<4;j++){						//updarting 4 conservative variable
					dw[cellNum][j]=-totalFlux[cellNum][j]*dt[cellNum]/area[j];
					w[cellNum][j]=wOld[cellNum][j]+alpha[rk]*dw[cellNum][j];
				}
				w[cellNum][4]=(gamma-1)*( w[cellNum][3] - 0.5*(w[cellNum][1]*w[cellNum][1]+w[cellNum][2]*w[cellNum][2])/w[cellNum][0]);						//updating pressure
			} //end of cell loop
		} //end of r-k loop

		switch (itr){						//prepare scaling of residues at the first iteration and print them
		case 1:	rp->FirstIteration();
				break;
		}

		if(itr%nfres==0 && itr!=1)			//calculate scaled residues and printing them
			rp->CalculateResidue(itr);

	} //end of iteration loop

	cout<<"\nIteration is finished with "<<itr<<" cycles."<<endl;;
	if(prntsln==1){
		cout<<"\nCell#"<<"     ro     "<<"      u      "<<"      v      "<<"     |V|     "<<"       P      "<<endl;
		cout<<"----------------------------------------------------------------------"<<endl;
		for(int i=0;i<celltot;i++)
		{
			cout<<setw(7)<<fixed<<left<<i<<setw(13)<<setprecision(3)<<fixed<<left<<w[i][0]<<setw(13)<<setprecision(3)<<fixed<<left<<w[i][1]/w[i][0]<< \
					setw(13)<<setprecision(3)<<fixed<<left<<w[i][2]/w[i][0]<<setw(13)<<setprecision(3)<<fixed<<left<<\
					sqrt(w[i][1]/w[i][0]*w[i][1]/w[i][0]+w[i][2]/w[i][0]*w[i][2]/w[i][0])<<setw(13)<<setprecision(3)<<fixed<<left<<w[i][4]<<endl;
		}
	}
}

double** Solver::GetW()
{
	return w;
}

