/*=================================================
               Cell Flux Class header
=================================================*/

#ifndef ONECELL_H_
#define ONECELL_H_

class Mesh;
class Utility;

class OneCell
{
public:
	OneCell(Mesh*,double**,double**,Utility*);
	void ClaculateCellFlux(int);
	void XYtoNT(double*,double*);
	void NTtoXY();
	void WallFlux();

private:
	Mesh* mp=nullptr;
	Utility* up=nullptr;
	double** w=nullptr;
	double** totalFlux=nullptr;
	int** connectivity=nullptr;
	int edgeNumber=0;
	double** nx=nullptr;
	double** ny=nullptr;
	double** dl=nullptr;
	double* wCellCurr=nullptr;
	double* wCellNeigh=nullptr;
	double* wFaceCurr=nullptr;
	double* wFaceNeigh=nullptr;
	double* faceFlux=nullptr;
	int cN=0;
	int fN=0;
	double fMomN=0;
	double fMomT=0;
	int neigh=0;
};

#endif
