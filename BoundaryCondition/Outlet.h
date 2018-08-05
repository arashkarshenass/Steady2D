/*=================================================
                Outlet Boundary Class header
=================================================*/
#ifndef OUTLET_H_
#define OUTLET_H_

class Outlet
{
public:
	void Initializer (double*,double*);
	void CalculateGhost();

private:
	double* wFaceCurr=nullptr;
	double* wGhost=nullptr;
};

#endif
