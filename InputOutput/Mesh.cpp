/*=================================================
                     mesh class source
=================================================*/

#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <string>

#include "Mesh.h"
#include "../Utility/Utility.h"
#include "../Inputs.h"
using namespace std;


void Mesh::Reader(char pub_name_in_func[])
{
	//-------opening mesh file----------
	for(int i=0;i<100;i++)
		fileName[i]=pub_name_in_func[i];

	FILE *pfile;							//Reading file
		pfile=fopen(fileName,"r");

		while(pfile==NULL)						//file is found or NOT?
		{
			printf("File cannot be found\n");
			printf("Enter file name: ");
			cin>>fileName;   		//or you can use cin>>filename;
			pfile=fopen(fileName,"r");
		}

	printf("\n%s is found. Reading mesh file ...\n",fileName);


	int type=0;

	//reading celltot
	fscanf(pfile,"%*72c %i",&celltot);

	//reading type
	fscanf(pfile,"%i",&type);

	//determining number of edges
	switch(type)
	{
		case 5: edgeNumber=3;
				break;
		case 9:	edgeNumber=4;
				break;
	}

	//allocating cellmap
	cellmap=new int*[celltot];
	for(int i=0;i<celltot;i++)
		cellmap[i]=new int[edgeNumber];

	//reading first cell
	switch(edgeNumber)
	{
	case 3:	fscanf(pfile,"%i %i %i %*i",*(cellmap+0)+0,*(cellmap+0)+1,*(cellmap+0)+2);
			break;
	case 4: fscanf(pfile,"%i %i %i %i %*i",*(cellmap+0)+0,*(cellmap+0)+1,*(cellmap+0)+2,*(cellmap+0)+3);
			break;
	}

	//reading rest of cells
	switch(edgeNumber)
	{
	case 3: for(int i=1;i<celltot;i++)
				fscanf(pfile,"%*i %i %i %i %*i",*(cellmap+i)+0,*(cellmap+i)+1,*(cellmap+i)+2);
			break;
	case 4: for(int i=1;i<celltot;i++)
				fscanf(pfile,"%*i %i %i %i %i %*i",*(cellmap+i)+0,*(cellmap+i)+1,*(cellmap+i)+2,*(cellmap+i)+3);
			break;
	}

	//reading node number
	fscanf(pfile,"%*31c %i",&nodetot);

	//allocating node cordinates
	x=new double[nodetot];
	y=new double[nodetot];

	//reading node cordinates
	for(int i=0;i<nodetot;i++)
		fscanf(pfile,"%lf %lf %*i",x+i,y+i);

	//reading number of boundary conditions
	fscanf(pfile,"%*31c %i",&bouTypes);

	//reading and allocating boundary nodes;
	bouPoints=new int**[bouTypes];
	bouSizes=new int[bouTypes];
	for(int i=0;i<bouTypes;i++)
	{
		int nn=0;
		fscanf(pfile,"%*s %*s");
		fscanf(pfile,"%*s %i",&nn);
		bouSizes[i]=nn;
		bouPoints[i]=new int*[nn];
		for(int j=0;j<nn;j++)
			bouPoints[i][j]=new int[2];
		for(int j=0;j<nn;j++)
			fscanf(pfile,"%*i %i %i",&bouPoints[i][j][0],&bouPoints[i][j][1]);
	}

	fclose(pfile);				//closing mesh file


	/*-------------------------print---------------------*/
	printf("Total cell number is %i      number of edges for each cell is %i      total node number is %i\n",celltot,edgeNumber,nodetot);
	printf("There exist %i type of boundary conditions in total.\n",bouTypes);


	if(prntmsh==1){
		puts("\n   Cell Mapping");
		puts("----------------------");
		for (int i=0;i<celltot;i++)			//printing cell nodes
		{	cout<<i<<"   ";
			for(int j=0;j<edgeNumber;j++)
				cout<<left<<setw(6)<<cellmap[i][j];
			cout<<endl;
		}

		puts("\n Nodes Coordinate");
		puts("----------------------");
		for (int i=0;i<nodetot;i++)			//printing node coordinates
			cout<<i<<"   "<<left<<setw(20)<<setprecision(16)<<fixed<<x[i]<<left<<setw(20)<<setprecision(16)<<fixed<<y[i]<<endl;
	}
}

void Mesh::Geometry(Utility* up)
{
	//allocating memory for dl,normal,xc,yc,xm,ym,area,minDis

	double dx[edgeNumber]={0};
	double dy[edgeNumber]={0};
	double sumx=0; double sumy=0;

	xc=new double[celltot];
	yc=new double[celltot];
	area=new double[celltot];
	minDis=new double[celltot];

	dl=new double*[celltot];
	for(int i=0;i<celltot;i++)
		dl[i]=new double[edgeNumber];

	normx=new double*[celltot];
	normy=new double*[celltot];
	for(int i=0;i<celltot;i++)
	{
		normx[i]=new double[edgeNumber];
		normy[i]=new double[edgeNumber];
	}


	for(int i=0;i<celltot;i++)
	{
		//calculating dl,nx and ny
		sumx=0;sumy=0;
		for(int j=0;j<edgeNumber;j++)
		{
			int j2=(j+1)%edgeNumber;
			dx[j]=x[cellmap[i][j2]]-x[cellmap[i][j]];		//delta x of edge j of cell i
			dy[j]=y[cellmap[i][j2]]-y[cellmap[i][j]];		//delta y of edge j of cell i
			dl[i][j]=sqrt(pow(dx[j],2)+pow(dy[j],2));				//length of edge j of cell i
			normx[i][j]=-dy[j]/dl[i][j];							//normal=nx i + ny j. here normal[0] matrix is nx
			normy[i][j]=dx[j]/dl[i][j];							//normal[1] matrix is ny
			sumx+=x[cellmap[i][j]];
			sumy+=y[cellmap[i][j]];
		}

		//calculating xc and yc
		xc[i]=sumx/edgeNumber;		//x coordinate of centroid of cell i
		yc[i]=sumy/edgeNumber;		//y coordinate of centroid of cell i

		//calculating area
		switch(edgeNumber)
		{
		case 3:
			area[i]=abs(0.5*up->VectorProduct2D(dx[0],dy[0],dx[1],dy[1]));
			break;
		case 4:
			area[i]=abs(0.5*up->VectorProduct2D(dx[0],dy[0],dx[1],dy[1]))+abs(0.5*up->VectorProduct2D(dx[2],dy[2],dx[3],dy[3]));
		}
	}

	//calculating minimum distance
	double rmin=1000;
	switch(edgeNumber)
	{
	case 3: for(int i=0;i<celltot;i++)
			{
				rmin=1000;
				for(int j=0;j<3;j++)
					if(dl[i][j]<rmin)
						rmin=dl[i][j];
				minDis[i]=0.7*rmin;
			}
			break;

	case 4:
			double rds[4]={0};
			for(int i=0;i<celltot;i++)
			{
				rmin=1000;
				rds[0]=sqrt(pow((x[cellmap[i][0]]+x[cellmap[i][1]])/2-xc[i],2)+pow((y[cellmap[i][0]]+y[cellmap[i][1]])/2-yc[i],2));
				rds[1]=sqrt(pow((x[cellmap[i][1]]+x[cellmap[i][2]])/2-xc[i],2)+pow((y[cellmap[i][1]]+y[cellmap[i][2]])/2-yc[i],2));
				rds[2]=sqrt(pow((x[cellmap[i][2]]+x[cellmap[i][3]])/2-xc[i],2)+pow((y[cellmap[i][2]]+y[cellmap[i][3]])/2-yc[i],2));
				rds[3]=sqrt(pow((x[cellmap[i][3]]+x[cellmap[i][0]])/2-xc[i],2)+pow((y[cellmap[i][3]]+y[cellmap[i][0]])/2-yc[i],2));
				for(int j=0;j<4;j++)
					if(rds[j]<rmin)
						rmin=rds[j];
				minDis[i]=rmin;
			}
			break;
	}

	double areaTotal=0;
	for(int i=0;i<celltot;i++)
		areaTotal+=area[i];

	//----------Print----------------

	if(prntmsh==1){
		puts("\n----delta L-------");
		up->Monitoring2D(celltot,edgeNumber,dl);

		puts("\n----normal X-------");
		up->Monitoring2D(celltot,edgeNumber,normx);

		puts("\n----normal Y-------");
		up->Monitoring2D(celltot,edgeNumber,normy);

		puts("\nCellCenterX    CellCenterY    Area     MinimumDistance");
		for(int i=0;i<celltot;i++)
			cout<<i<<"  "<<xc[i]<<"  "<<yc[i]<<"  "<<area[i]<<"  "<<minDis[i]<<endl;
		cout<<"Total area is: "<<areaTotal<<endl;
	}
}


void Mesh::Connectivity()
{
	//allocating memory for connectivity matrix
	neigh=new int*[celltot];
	for(int i=0;i<celltot;i++)
		neigh[i]=new int[edgeNumber];

	//initializing connectivity matrix
	for(int i=0;i<celltot;i++)
		for(int j=0;j<edgeNumber;j++)
			neigh[i][j]=-1;

	int flag,m,n;
	for(int i=0;i<celltot;i++)
	{
		for(int j=0;j<edgeNumber;j++)
		{
			if(neigh[i][j]!=-1)
				continue;
			flag=0;
			for(int k=i+1;k<celltot;k++)
			{
				if(flag==1)
					break;

				for(int l=0;l<edgeNumber;l++)
				{
					if(l==j && k==i)
						continue;

					if(cellmap[i][j]==cellmap[k][l])
					{
						switch(edgeNumber)
						{
						case 4:	if(j==3)
									m=0;
								else
									m=j+1;

								if(l==0)
									n=3;
								else
									n=l-1;
								break;

						case 3: if(j==2)
									m=0;
								else
									m=j+1;

								if(l==0)
									n=2;
								else
									n=l-1;
								break;
						}

						if(cellmap[i][m]==cellmap[k][n])
						{
							flag=1;
							neigh[i][j]=k;
							neigh[k][n]=i;
							break;
						}
					}
				}
			}
		}
	}

	//find B.C surfaces from bouPoints
	for(int i=0;i<2;i++)		//it is set to 2 because zero is for inlet and 1 is for outlet. unchaged neigh is wall by default
	{
		for(int j=0;j<bouSizes[i];j++)
		{
			int a=bouPoints[i][j][0];int b=bouPoints[i][j][1];
			for(int k=0;k<celltot;k++)
			{
				for(int l=0;l<edgeNumber;l++)
				{
					if(cellmap[k][l]==a)
					{
						if(l==edgeNumber-1)
						{
							if(cellmap[k][0]==b)
							{
								switch(i)
								{
								case 0:neigh[k][l]=-2;
									   break;
								case 1:neigh[k][l]=-3;
									   break;
								}
							}
							else if(cellmap[k][l-1]==b)
							{
								switch(i)
								{
								case 0:neigh[k][l-1]=-2;
										break;
								case 1:neigh[k][l-1]=-3;
										break;
								}
							}
						}

						else if(l==0)
						{
							if(cellmap[k][1]==b)
							{
								switch(i)
								{
								case 0:neigh[k][l]=-2;
										break;
								case 1:neigh[k][l]=-3;
										break;
								}
							}
							else if(cellmap[k][edgeNumber-1]==b)
							{
								switch(i)
								{
								case 0:neigh[k][edgeNumber-1]=-2;
										break;
								case 1:neigh[k][edgeNumber-1]=-3;
										break;
								}
							}
						}

						else
						{
							if(cellmap[k][l+1]==b)
							{
								switch(i)
								{
								case 0:neigh[k][l]=-2;
									   break;
								case 1:neigh[k][l]=-3;
										break;
								}
							}
							else if(cellmap[k][l-1]==b)
							{
								switch(i)
								{
								case 0:neigh[k][l-1]=-2;
									break;
								case 1:neigh[k][l-1]=-3;
									break;
								}
							}

						}
					}
				}

			}
		}
	}


	if(prntmsh==1){
		puts("\n    Connectivity");
		puts("----------------------");
		for(int i=0;i<celltot;i++)
			cout<<i<<"  "<<left<<setw(6)<<neigh[i][0]<<left<<setw(6)<<neigh[i][1]<<left<<setw(6)<<neigh[i][2]<<left<<setw(6)<<neigh[i][3]<<endl;

		puts("\nConnectivity is created successfully");
	}
}


int Mesh::GetCelltot() const
{
	return celltot;
}

int Mesh::GetNodetot() const
{
	return nodetot;
}


int Mesh::GetEdgeNumber() const
{
	return edgeNumber;
}


int** Mesh::GetCellmap() const
{
	return cellmap;
}


double* Mesh::GetX() const
{
	return x;
}


double* Mesh::GetY() const
{
	return y;
}


double** Mesh::GetDl() const
{
	return dl;
}


double** Mesh::GetNx() const
{
	return normx;
}

double** Mesh::GetNy() const
{
	return normy;
}


double* Mesh::GetArea() const
{
	return area;
}


double* Mesh::GetMinDis() const
{
	return minDis;
}


int** Mesh::GetConnectivity() const
{
	return neigh;
}

