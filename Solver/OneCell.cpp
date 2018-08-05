/*=================================================
                Flux Class source
=================================================*/


#include "OneCell.h"

#include "../InputOutput/Mesh.h"
#include "../Inputs.h"
#include "../Utility/Utility.h"
#include "../BoundaryCondition/Inlet.h"
#include "../BoundaryCondition/Outlet.h"
#include "../BoundaryCondition/Wall.h"
#include "../FaceFlux/StegerWarmingVector.h"
#include "../FaceFlux/StegerWarmingmatrix.h"
#include "../FaceFlux/Upwind.h"
#include "../FaceFlux/AUSM.h"
#include "../FaceFlux/Average.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdlib.h>

using namespace std;

Inlet inlet;													//create inlet object
Outlet outlet;													//create outlet object
Wall wall;														//create wall object
AUSM ausm;														//create ausm object
StegerWarmingVector swv;
StegerWarmingMatrix swm;
Upwind upwind;
Average average;

OneCell::OneCell(Mesh* mpTemp,double** wTemp,double** totalFluxTemp,Utility* upTemp){
	mp=mpTemp;													//mesh object pointer from solver obj
	w=wTemp;													//conservative variables of all cell from solver obj
	totalFlux=totalFluxTemp;									//total flux of all cells coming from solver obj
	up=upTemp;													//utility object pointer from solver

	connectivity=mp->GetConnectivity();							//get connectuvity of cells from mesh obj
	edgeNumber=mp->GetEdgeNumber();								//get edge number of cells frommesh obj
	nx=mp->GetNx();												//get normal x from mesh object
	ny=mp->GetNy();												//get normal y from mesh obj
	dl=mp->GetDl();												//get face length from mesh obj

	wCellCurr=new double[5];									//current cell conservative variables in x-y coord.
	wCellNeigh=new double[5];									//neighbr cell conservative variables in x-y coord.
	wFaceCurr=new double[5];									//variables of face from current cell side in n-t coord.
	wFaceNeigh=new double[5];									//variables of face from neighbr cell side in n-t coord.
	faceFlux=new double[4];										//fluxes passing a face

	inlet.Initializer(wFaceCurr,wCellNeigh);					//initializing inlet object
	outlet.Initializer(wFaceCurr,wCellNeigh);					//initializing outletaksel object
	wall.Initializer(wFaceCurr,wFaceNeigh);						//initializing wall object
	swv.Initializer(wFaceCurr,wFaceNeigh,faceFlux,dl);			//initializing steger object
	swm.Initializer(wFaceCurr,wFaceNeigh,faceFlux,dl);			//initializing steger object
	upwind.Initializer(wFaceCurr,wFaceNeigh,faceFlux,dl);		//initializing upwind object
	ausm.Initializer(wFaceCurr,wFaceNeigh,faceFlux,dl);			//initializing ausm object
	average.Initializer(wFaceCurr,wFaceNeigh,faceFlux,dl);		//initializing average object
}

void OneCell::ClaculateCellFlux(int cellNum){
	cN=cellNum;
	wCellCurr[0]=w[cN][0];wCellCurr[1]=w[cN][1];			//variables of current cell (in x-y coord) is w[n]
	wCellCurr[2]=w[cN][2];wCellCurr[3]=w[cN][3];
	wCellCurr[4]=w[cN][4];


	for(fN=0;fN<edgeNumber;fN++){						//loop for faces of a cell
		XYtoNT(wCellCurr,wFaceCurr);					//convert variables of current cell from xy to nt coordinate
		neigh=connectivity[cN][fN];						//neighbour cell number

		switch(neigh){
		case -1:
				wall.CalculateGhost();						//calculates face variable in n-t from wall ghost cell side and store in
				break;										//wFaceNeigh when wall boundary exist
		case -2:
				inlet.CalculateGhost();					//calculates variables at inlet ghost cell in x-y and store in wCellNeigh
				XYtoNT(wCellNeigh,wFaceNeigh);				//convert variables of inlet cell from xy to nt coordinate
				break;
		case -3:
				outlet.CalculateGhost();					//calculates variables at outlet ghost cell in x-y and store in wCellNeigh
				XYtoNT(wCellNeigh,wFaceNeigh);				//convert variables of outlet cell from xy to nt coordinate
					break;
		default:
				wCellNeigh[0]=w[neigh][0];					//variables at interior neighbor cell in x-y and store in wCellNeigh
				wCellNeigh[1]=w[neigh][1];
				wCellNeigh[2]=w[neigh][2];
				wCellNeigh[3]=w[neigh][3];
				wCellNeigh[4]=w[neigh][4];
				XYtoNT(wCellNeigh,wFaceNeigh);				//convert variables of neighbor cell from xy to nt coordinate
				break;
		}	// end of switch												//end of switch

		switch(advectionType){
		case 1:	//averaging
				average.CalculateFaceFlux(cN,fN);			//calculate face flux at nt using average in vector form
				break;
		case 2:	//upwind
				upwind.CalculateFaceFlux(cN,fN);			//calcualte face flux at nt using upwind in vector form
				break;
		case 3:	//sw vector
				swv.CalculateFaceFlux(cN,fN);		//calculate face flux by stegerwarming in vector form
				break;
		case 4:	//sw matrix
				swm.CalculateFaceFlux(cN,fN);				//calculate face flux by steger warming in matrix form
				break;
		case 5:	//ausm
				ausm.CalculateFaceFlux(cN,fN);			//calculate face flux at nt using ausm	in vector form
				break;
		default:
				cout<<"Error in advection type"<<endl;
				exit (EXIT_FAILURE);
				break;
		}

		NTtoXY();												//nt to xy back transfer
		totalFlux[cN][0]+=faceFlux[0];
		totalFlux[cN][1]+=faceFlux[1];
		totalFlux[cN][2]+=faceFlux[2];
		totalFlux[cN][3]+=faceFlux[3];
	}	//end of face loop
}	//end of function


void OneCell::XYtoNT(double* wCell,double* wFace){
	wFace[0]=wCell[0];
	wFace[1]=wCell[1]*nx[cN][fN]+wCell[2]*ny[cN][fN];
	wFace[2]=-wCell[1]*ny[cN][fN]+wCell[2]*nx[cN][fN];
	wFace[3]=wCell[3];
	wFace[4]=wCell[4];
}

void OneCell::NTtoXY(){
	fMomN=faceFlux[1]; fMomT=faceFlux[2];
	faceFlux[1]=fMomN*nx[cN][fN]-fMomT*ny[cN][fN];
	faceFlux[2]=fMomN*ny[cN][fN]+fMomT*nx[cN][fN];
}

void OneCell::WallFlux(){
	faceFlux[0]=0;
	faceFlux[1]=wFaceCurr[4]*dl[cN][cN];
	faceFlux[2]=0;
	faceFlux[3]=0;
}
