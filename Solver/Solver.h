/*=================================================
                Solver Class header
=================================================*/

#ifndef SOLVER_H_
#define SOLVER_H_

class Mesh;
class LocalTimeStep;
class OneCell;
class Residue;
class Utility;


//////////////////////////////
class Solver {
public:
	Solver(Mesh*,Utility*);
	void Iterator();
	double** GetW();

private:
	int celltot=0;						//total cell number
	double **w=nullptr;					//pointer to conservative variables, w
	double **wOld=nullptr;				//pointer to wOld
	double **dw=nullptr;				//pointer to change of cons. variables in a cell, dw, in heap
	double **totalFlux=nullptr;			//pointer to total fluxes
	double *dt=nullptr;					//pointer to local time steps
	double *area=nullptr;				//pointer to cell area
	double* res=nullptr;				//pointer to residues
	int		itr=0;						//iteration counter
	Mesh* mp=nullptr;					//pointer to mesh object
	LocalTimeStep* tp=nullptr;
	OneCell* ocp=nullptr;
	Residue* rp=nullptr;
	Utility* up=nullptr;
};

#endif
