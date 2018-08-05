/*=================================================
                Residue Class source
=================================================*/
#include "Residue.h"
#include "../Utility/Utility.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
using namespace std;

Residue::Residue(double** dwTemp,double* resTemp, Utility* upTemp,int celltotTemp){
	celltot=celltotTemp;
	dw=dwTemp;
	res=resTemp;
	dwFirst=new double[4];
	up=upTemp;


}

void Residue::FirstIteration(){

	up->VectorInitializer(res,4,0);				//set residues to 0
	double tot[4]={0};
	for (int i=0;i<4;i++){						//average resudues of the 1st iteration to be used for scaling of next iterations
		for (int j=0;j<celltot;j++)
			tot[i]=tot[i]+dw[j][i]*dw[j][i];
		dwFirst[i]=sqrt(tot[i])/celltot;
		res[i]=log10(dwFirst[i]);
	}

	cout<<setw(10)<<fixed<<left<<1<<setw(13)<<setprecision(8)<<fixed<<left<<res[0]<<setw(13)<<setprecision(8)<<fixed<<left<<res[1]<< \
			setw(13)<<setprecision(8)<<fixed<<left<<res[2]<<setw(13)<<setprecision(8)<<fixed<<left<<res[3]<<endl;
}

void Residue::CalculateResidue(int itr){

	for (int i=0;i<celltot;i++)							//convert residues to scaled residue
		for (int j=0;j<4;j++)
			dw[i][j]=abs(dw[i][j])/dwFirst[j];

	up->VectorInitializer(res,4,0);						//set residues to 0
	for(int i=0;i<4;i++){								//find maximum scaled residue
		for(int j=0;j<celltot;j++)
			if(dw[j][i]>res[i])
				res[i]=dw[j][i];
		res[i]=log10(res[i]);
	}
	cout<<setw(10)<<fixed<<left<<itr<<setw(13)<<setprecision(8)<<fixed<<left<<res[0]<<setw(13)<<setprecision(8)<<fixed<<left<<res[1]<< \
		setw(13)<<setprecision(8)<<fixed<<left<<res[2]<<setw(13)<<setprecision(8)<<fixed<<left<<res[3]<<endl;
}
