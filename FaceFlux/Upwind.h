/*=================================================
                Steger Class header
=================================================*/
#ifndef UPWIND_H_
#define UPWIND_H_

class Upwind {
public:
	void Initializer(double*,double*,double*,double**);
	void CalculateFaceFlux(int,int);


private:
	double* wFaceCurr=nullptr;
	double* wFaceNeigh=nullptr;
	double* faceFlux=nullptr;
	double** dl=nullptr;
};

#endif
