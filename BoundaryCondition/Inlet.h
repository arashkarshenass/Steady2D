/*=================================================
                Inlet Boundary Class header
=================================================*/
#ifndef INLET_H_
#define INLET_H_

class Inlet
{
public:
	void Initializer (double*,double*);
	void CalculateGhost();

private:
	double* wFaceCurr=nullptr;
	double* wCellNeigh=nullptr;
};

#endif
