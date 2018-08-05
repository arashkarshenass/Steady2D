/*=================================================
                     LocalTimeStep class header
=================================================*/

#ifndef LOCALTIMESTEP_H_
#define LOCALTIMESTEP_H_
class Mesh;

class LocalTimeStep {
public:
	LocalTimeStep(double**,Mesh*,double*,int);
	void ClaculateLTS();

private:
	double* dt=nullptr;
	double** w=nullptr;
	double* minDis=nullptr;
	int celltot=0;
	double u=0;
	double v=0;
	double flowSpeed=0;
	double soundSpeed=0;
	double totalSpeed=0;

};



#endif
