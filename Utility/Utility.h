

#ifndef UTILITY_UTILITY_H_
#define UTILITY_UTILITY_H_

class Utility {
public:
	void VectorInitializer(double*,int,double);
	void MatrixInitializer(double**,int,int,double);
	double VectorProduct2D(double,double,double,double);
	void VectorEqualizer(double*,double*,int);
	void MatrixEqualizer(double**,double**,int,int);
	void Monitoring2D(int,int,double**);
private:
	double m=0;
	int colN=0;
	int size=0;
};

#endif
