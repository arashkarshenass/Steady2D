#include <iostream>
#include <iomanip>
#include "Utility.h"

using namespace std;

void Utility::VectorInitializer(double* vector, int size, double value)
{
	for(int i=0;i<size;i++)
		vector[i]=value;
}

void Utility::MatrixInitializer(double** matrix,int rowN,int colN,double value)
{
	for(int i=0;i<rowN;i++)
	{
		for(int j=0;j<colN;j++)
		{
			matrix[i][j]=value;
		}
	}
}

double Utility::VectorProduct2D(double x1, double y1, double x2, double y2)
{
	return m=x1*y2-y1*x2;
}

void Utility::VectorEqualizer(double* target, double* ref,int size)
{
	for(int i=0;i<size;i++)
		target[i]=ref[i];
}

void Utility::MatrixEqualizer(double** target, double** ref, int rowN,int colN)
{
	for(int i=0;i<rowN;i++)
		for(int j=0;j<colN;j++)
			target[i][j]=ref[i][j];
}


void Utility::Monitoring2D(int cellTot, int edgeNumber, double** matrix)
{
	for(int i=0;i<cellTot;i++)
		{	cout<<i<<"  ";
			for(int j=0;j<edgeNumber;j++)
				cout<<setw(17)<<setprecision(15)<<fixed<<left<<matrix[i][j]<<"  ";
			cout<<endl;
		}
}

