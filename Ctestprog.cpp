#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <mpi.h>

using namespace std;


int main(int argc, char **argv)

{
	const int nx=250, ny=250, nz=30;
	const int ndim=(nx*ny*nz);
	const float xsize=1000, ysize=1000, zsize=150;
	int numXNodes, numUNodes, numElements, numLocalNodes;

	int density;

	float *xNode, *yNode, *zNode;
	float *uNode, *vNode, *wNode;

	float *xAbsorb, *yAbsorb, *zAbsorb;
	float *xSource, *ySource, *zSource;

	float t, dt, ltime, x,y,z, dx, dy, dz;
	int i, j, k, n;


	int currentNonLinearIteration, numIterations;
	bool remesh;

	int ierr, mpiID;
	// char logfile[256], mpich;
	string logname, mpich;


	int ncpus;
	float offsetx, largedx;

	xNode=new float[ndim];
	yNode=new float[ndim];
	zNode=new float[ndim];

	uNode=new float[ndim];
	vNode=new float[ndim];
	wNode=new float[ndim];

	xAbsorb=new float[ndim];
	yAbsorb=new float[ndim];
	zAbsorb=new float[ndim];

	xSource=new float[ndim];
	ySource=new float[ndim];
	zSource=new float[ndim];	

	cout << "begin testprog\n";

//	void MPI_Init[ierr];
//	void MPI_Comm_size[MPI_COMM_WORLD][ncpus][ierr];
//	void MPI_Comm_rank[MPI_COMM_WORLD][mpiID][ierr];

	ierr=MPI_Init(&argc, &argv);
	ierr=MPI_Comm_size(MPI_COMM_WORLD, &ncpus);
	ierr=MPI_Comm_rank(MPI_COMM_WORLD, &mpiID);
	
	largedx = 1000./ncpus;
	offsetx = largedx * mpiID;

	stringstream logstr;

	logstr << "output-"<< mpiID << ".log";
	logname = logstr.str();

	//ifstream OpenFile(678, file=logfile, status='replace')
	
	std::ofstream logfile;
	
	logfile.open(logname.c_str());
	
	
	numXNodes=ndim;
	numUNodes=ndim;

	dx=xsize/(nx-1);
	dy=ysize/(ny-1);
	dz=zsize/(nz-1);


	n=1;
	for ( int k=0; k<= nz-1; k=k+1)
	{
		for ( int j=0; j<=ny-1; j=j+1)
		{
			for ( int i=0; i<=nx-1; i=i+1)
			{
				x= i*dx + offsetx;
				y= j*dy;
				z= k*dz;

				xNode[n]=x;
				yNode[n]=y;
				zNode[n]=z;
				uNode[n]=5.;
				vNode[n]=0.;
				wNode[n]=0.;

				n=n+1;
			}
		}
	}
	
	density=1.227;
  
	t=0;
	dt=1;
	remesh=false;
	currentNonLinearIteration=1;
	numIterations=2;

	ltime = 100000;

	for (int i=1; i<=int(ltime/dt)+1; i=i+1)
	{
		for (int j=1; j<2; j=j+1)
		{
			void turbineFarmFluidityInterface(
             int numXNodes, int numUNodes,
             float *xNode, float *yNode, float *zNode, 
             float *uNode, float *vNode, float *wNode, 
             float density, 
             float *xAbsorb, float *yAbsorb, float *zAbsorb, 
             float *xSource, float *ySource, float *zSource, 
             float t, float dt, 
             int currentNonLinearIteration, int numIterations, 
             bool remesh );
		}
     t=t+dt;
	}
	
	MPI_Finalize();

	delete [] xNode;
	delete [] yNode;
	delete [] zNode;

	delete [] uNode;
	delete [] vNode;
	delete [] wNode;

	delete [] xAbsorb;
	delete [] yAbsorb;
	delete [] zAbsorb;

	delete [] xSource;
	delete [] ySource;
	delete [] zSource;
  
	cout << "end testprog\n";
	
	return 0;	
}
