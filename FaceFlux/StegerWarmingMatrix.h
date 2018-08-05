
#ifndef FACEFLUX_STEGERWARMINGMATRIX_H_
#define FACEFLUX_STEGERWARMINGMATRIX_H_

class StegerWarmingMatrix {
public:
	void Initializer(double*,double*,double*,double**);
	void CalculateFaceFlux(int,int);
	void CalculateacobianMatrix();

private:
	double* wFaceCurr=nullptr;
	double* wFaceNeigh=nullptr;
	double* faceFlux=nullptr;
	double** dL=nullptr;
	double g1=0;
	double* lam=nullptr;
	double un=0;
	double ut=0;
	double a=0;
	double s2=0;
	double** jacob=0;
};

#endif
