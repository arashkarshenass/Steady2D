/*=================================================
                Average Class header
=================================================*/
#ifndef FACEFLUX_AVERAGE_H_
#define FACEFLUX_AVERAGE_H_

class Average
{public:
	void Initializer(double*,double*,double*,double**);
	void CalculateFaceFlux(int,int);

private:
	double* wFaceCurr=nullptr;
	double* wFaceNeigh=nullptr;
	double* faceFlux=nullptr;
	double** dL=nullptr;
};

#endif
