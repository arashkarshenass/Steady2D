/*=================================================
                Wall Boundary Class header
=================================================*/
#ifndef WALL_H_
#define WALL_H_

class Wall
{
public:
	void Initializer(double*,double*);
	void CalculateGhost();

private:
	double* wFaceCurr=nullptr;
	double* wFaceNeigh=nullptr;
};

#endif
