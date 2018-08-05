/*=================================================
                AUSM Class header
=================================================*/
#ifndef FACEFLUX_AUSM_H_
#define FACEFLUX_AUSM_H_

class AUSM
{
public:
	void Initializer(double*,double*,double*,double**);
	void CalculateFaceFlux(int,int);
private:
	double* wFaceCurr=nullptr;
	double* wFaceNeigh=nullptr;
	double* faceFlux=nullptr;
	double** dL=nullptr;
	double* lam=nullptr;
	double* fluxP=nullptr;
	double* fluxM=nullptr;
	double  ro=0;
	double  u=0;
	double  v=0;
	double  p=0;
	double  a=0;
	double  m=0;
	double  mP=0;
	double  mM=0;
	double  pP=0;
	double  pM=0;
};

#endif /* FACEFLUX_AUSM_H_ */
