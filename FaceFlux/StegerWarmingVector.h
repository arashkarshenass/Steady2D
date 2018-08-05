/*=================================================
                Steger Class header
=================================================*/
#ifndef STEGERWARMINGVECTOR_H_
#define STEGERWARMING_H_

class StegerWarmingVector {
public:
	void Initializer(double*,double*,double*,double**);
	void CalculateFaceFlux(int,int);

private:
	double* wFaceCurr=nullptr;
	double* wFaceNeigh=nullptr;
	double* faceFlux=nullptr;
	double** dL=nullptr;
	double 	g1=0;
};

#endif
