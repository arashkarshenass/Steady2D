/*
 ==================================================
                2D Steady Euler in NT coordinate
 ==================================================
   1D Euler Flow Solver on unstructured Grids
                  developed by
               Arash Karshenass
       Middle East Technical University
                 Ankara- TURKEY
 ================================================== */

#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <ctime>
#include "InputOutput/Mesh.h"
#include "InputOutput/OutPut.h"
#include "Solver/Solver.h"
#include "Utility/Utility.h"
#include "Inputs.h"

using namespace std;

//--------Inputs----------
//through input header file


int main()
{
	clock_t time_start=clock();
	printf("------Program is started-----\n");
	char fileName[]="NASA_2D_C++.su2";     						//name of mesh file

	Utility utility;
	Mesh mesh;
	mesh.Reader(fileName);
	mesh.Geometry(&utility);
	mesh.Connectivity();

	Solver solver(&mesh,&utility);

	clock_t time_end=clock();

	if(savesln==1){
	OutPut output(&mesh,&solver);
	output.TecPlot2D();
	}
	/*else{
		puts("You have not saved your results. Do you want to save them in Tecplot formast? Press ''y''  for YES and ''n'' for NO." );
		char save;
		cin>>save;
		int cont=0;
		while(cont==0){
			if(save=='y' || save=='Y'){
			OutPut output(&mesh,&solver);
			output.TecPlot2D();
			cont=1;
			}
			else if(save=='n' || save=='N')
				cont=1;
			else{
				puts("Wrong input. ''y'' or ''n''?");
				cin>>save;
			}
		}
	}*/

	puts("\nProgram has terminated successfully :)");
	cout<<"Run time: "<<(time_end - time_start)/(double)CLOCKS_PER_SEC<<endl;

	return 0;
}
