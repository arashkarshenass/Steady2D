
#ifndef OUTPUT_H_
#define OUTPUT_H_

class Mesh;
class Solver;

class OutPut{
public:
	OutPut(Mesh*,Solver*);
	void TecPlot2D();
private:
	double* x=nullptr;
	double* y=nullptr;
	int** cellmap=nullptr;
	int celltot=0;
	int nodetot=0;
	double** w=nullptr;
	double* rho=nullptr;
	double* u=nullptr;
	double* v=nullptr;
	double* p=nullptr;
	double* pTot=nullptr;
	double* t=nullptr;
	double* tTot=nullptr;
	double* m=nullptr;
	double* e=nullptr;
	double* eTot=nullptr;
	double* h=nullptr;
	double* hTot=nullptr;
	double* ds=nullptr;
};

#endif
