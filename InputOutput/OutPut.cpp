#include "OutPut.h"
#include "Mesh.h"
#include "../Inputs.h"
#include "../Solver/Solver.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

OutPut::OutPut(Mesh* mp,Solver* sp){
	celltot=mp->GetCelltot();
	nodetot=mp->GetNodetot();
	x=mp->GetX();
	y=mp->GetY();
	cellmap=mp->GetCellmap();
	w=sp->GetW();

	rho=new double[celltot];
	u=new double[celltot];
	v=new double[celltot];
	p=new double[celltot];
	pTot=new double[celltot];
	t=new double[celltot];
	tTot=new double[celltot];
	m=new double[celltot];
	e=new double[celltot];
	eTot=new double[celltot];
	h=new double[celltot];
	hTot=new double[celltot];
	ds=new double[celltot];
	double a;

	//calculating primitive variables
	for(int i=0;i<celltot;i++){
		rho[i]=w[i][0];												//density
		u[i]=w[i][1]/w[i][0];										//x-velocity
		v[i]=w[i][2]/w[i][0];										//y-velocity
		p[i]=w[i][4];												//pressure
		t[i]=p[i]/rho[i]/R;											//temperature
		e[i]=cv*t[i];												//internal energy
		eTot[i]=w[i][3]/w[i][0];									//total energy
		h[i]=cp*t[i];												//enthalpy
		a=sqrt(gamma*p[i]/rho[i]);									//speed of sound
		m[i]=sqrt(u[i]*u[i]+v[i]*v[i])/a;							//Mach number
		tTot[i]=t[i]*(1+(gamma-1)/2*m[i]*m[i]);						//Total temperature
		pTot[i]=p[i]*pow(1+(gamma-1)/2*m[i]*m[i],gamma/(gamma-1));	//Total pressure
		hTot[i]=cp*tTot[i];											//Total enthalpy
	}

	double tRef=t[0];
	double rhoRef=rho[0];
	for(int i=0;i<celltot;i++)
		ds[i]=cv*log(t[i]/tRef)+R*log(rhoRef/rho[i]);
}

void OutPut::TecPlot2D(){

	puts("\nSaving primitive cell centered solution in output file with TecPlot Format ...");

	ofstream slnfile;
	slnfile.open("solution.dat");
	slnfile<<"TITLE=\"Solution\"\nFILETYPE=FULL"<<endl;
	slnfile<<"VARIABLES=\"X\" \"Y\" \"Rho\" \"u\" \"v\" \"P\" \"P tot\" \"T\" \"T tot\" \"M\" \"e\" \"e tot\" \"h\" \"h tot\" \"dS\""<<endl;
	slnfile<<"ZONE\nT=Interior\nNODES="<<nodetot<<endl;
	slnfile<<"ELEMENTS="<<celltot<<endl;
	slnfile<<"DATAPACKING=BLOCK\nZONETYPE=FEQUADRILATERAL"<<endl;
	slnfile<<"VARLOCATION=([3-15]=CELLCENTERED)"<<endl;
	slnfile<<setprecision(6)<<fixed<<scientific;

	slnfile<<"#x-coord of vertices"<<endl;
	for(int i=0;i<nodetot;i++)
		slnfile<<x[i]<<endl;

	slnfile<<"#y-coord of vertices"<<endl;
	for(int i=0;i<nodetot;i++)
		slnfile<<y[i]<<endl;

	slnfile<<"#Rho of cell centers"<<endl;
	for(int i=0;i<celltot;i++)
		slnfile<<rho[i]<<endl;

	slnfile<<"#X-Velocity of cell centers"<<endl;
	for(int i=0;i<celltot;i++)
		slnfile<<u[i]<<endl;

	slnfile<<"#Y-Velocity of cell centers"<<endl;
	for(int i=0;i<celltot;i++)
		slnfile<<v[i]<<endl;

	slnfile<<"#Pressure of cell centers"<<endl;
	for(int i=0;i<celltot;i++)
		slnfile<<p[i]<<endl;

	slnfile<<"#Pressure Total of cell centers"<<endl;
	for(int i=0;i<celltot;i++)
		slnfile<<pTot[i]<<endl;

	slnfile<<"#Temperature of cell centers"<<endl;
	for(int i=0;i<celltot;i++)
		slnfile<<t[i]<<endl;

	slnfile<<"#Temperature total of cell centers"<<endl;
	for(int i=0;i<celltot;i++)
		slnfile<<tTot[i]<<endl;

	slnfile<<"#Mach of cell centers"<<endl;
	for(int i=0;i<celltot;i++)
		slnfile<<m[i]<<endl;

	slnfile<<"#Energy of cell centers"<<endl;
	for(int i=0;i<celltot;i++)
		slnfile<<e[i]<<endl;

	slnfile<<"#Energy Total of cell centers"<<endl;
	for(int i=0;i<celltot;i++)
		slnfile<<eTot[i]<<endl;

	slnfile<<"#Enthalpy of cell centers"<<endl;
	for(int i=0;i<celltot;i++)
		slnfile<<h[i]<<endl;

	slnfile<<"#Enthalpy Total of cell centers"<<endl;
	for(int i=0;i<celltot;i++)
		slnfile<<hTot[i]<<endl;

	slnfile<<setprecision(6)<<fixed;
	slnfile<<"#Entropy Change of cell centers"<<endl;
	for(int i=0;i<celltot;i++)
		slnfile<<ds[i]<<endl;

	slnfile<<"#Cellmap of cells"<<endl;
	for(int i=0;i<celltot;i++)
		slnfile<<cellmap[i][0]+1<<" "<<cellmap[i][1]+1<<" "<<cellmap[i][2]+1<<" "<<cellmap[i][3]+1<<endl;

	slnfile.close();
	puts("Data are saved successfully :)");
}
