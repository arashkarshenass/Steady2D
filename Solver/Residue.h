/*=================================================
                Residue Class header
=================================================*/

#ifndef RESIDUE_H_
#define RESIDUE_H_

class Solver;
class Utility;

class Residue
{
public:
	Residue(double**,double*, Utility*,int);
	void FirstIteration();
	void CalculateResidue(int);

private:
	int celltot=0;
	double** dw=nullptr;
	double* res=nullptr;
	double* dwFirst=nullptr;
	Utility* up=nullptr;
};

#endif
