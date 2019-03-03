/* Function Files for task 5
 *	
 * Written by: Ion Berasaluce
 * CID 00950097
 */

#include <mpi.h>
#include <iostream>
#include <fstream>

#include "headers.h"
using namespace std;

#define F77NAME(x) x##_
extern "C" {
	void F77NAME(dcopy)(const int& N, double* DX, const int& INCX, double* DY,
	    				const int& INCY);

	void F77NAME(pdpbsv)(const char* UPLO, const int* N, const int* k, const int* NRHS,  double* A, const int* ja, const int* desca,
			double* B, const int* ib, const int* descb, double* work, const int* lwork, int* INFO);
	
	void Cblacs_get(int, int, int*);
	void Cblacs_pinfo(int*, int*);
	void Cblacs_gridinit(int*, char*, int, int);
	void Cblacs_gridinfo(int, int*, int*, int*, int*);
	void Cblacs_exit(int);
	void Cblacs_gridexit(int);
	
}

void implicitSolvPar(const int subdomain_size, const int subdomain_start, double* K, double* M,
					 double dt, double T, double T1, double beta, const int DoF, double l){
/* subdomain size is the domain of the partitioned DoF of the specific process
   subdomain starts is the starting point of the Degrees of Freedom of the process
   DoF is the total number of Degrees of Freedom of the system excluding those at the boundary conditions
   M is the mass matrix of the system
   K is the assembled full stiffness matrix
   dt is the timestep used 
   T is the maximum time 
   T1 is the maximum loading time 
   l is the length of the elements
   beta is the newmark constant
*/
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	// Initialize the required variables for the solver
	double k1 = 1.0/(beta*dt*dt);         double k2 = 1.0/(dt*beta);              double k3 = (1.0/(2*beta) - 1); 
	double n[subdomain_size] = {0};       double ndot[subdomain_size] = {0};      double lambda = 0.5;
  	double nddot1[subdomain_size] = {0};  double nddot[subdomain_size] = {0};     double t = 0;
  	double F[subdomain_size];             double Fy,qy;  					      double Keff[6*subdomain_size];
  	double k4[DoF] = {0}; 

  	// Lapack solver requires the size of all of the procesess to be the same so we make a new variable 
  	// just for the solver to think they are the same.
  	int subsize_lappack = (double)(DoF+3)/(double)size;

  	// Set up ScaLapack solver stuff
  	char U = 'U';
	int NRHS = 1;
	int ja = 1;
	int ib = 1;
	int lda = 6;
	int k = 5;
	int ctx;
	int INFO = 0;
	// We know k is larger than the right hand side always
	int lwork = (subsize_lappack + 2*k)*k + k*k;  

	// Set up the CBLACS stuff
	char order = 'R';
	int nrow = 1;
	int ncol = size;
	int mype;
	int npe;
	int myrow;
	int mycol;
	Cblacs_pinfo(&mype , &npe);
	Cblacs_get(0, 0, &ctx);
	Cblacs_gridinit(&ctx, &order, 1, npe);
	Cblacs_gridinfo(ctx, &nrow, &ncol, &myrow, &mycol);

	// Descriptor array of the first array 
	int desca[7];
		desca[0] = 501;				// Type of banded Matrix 1-P
		desca[1] = ctx;				// Context
		desca[2] = DoF;				// Problem Size
		desca[3] = subsize_lappack;	// Blocking 
		desca[4] = 0;				// Process/row/column
		desca[5] = lda;				// Local Size
		desca[6] = 0;				// Reserved

	// Descriptor array of the first array 
	int descb[7];
		descb[0] = 502;				// Type of banded Matrix 1-P
		descb[1] = ctx;				// Context
		descb[2] = DoF;				// Problem Size
		descb[3] = subsize_lappack;	// Blocking 
		descb[4] = 0;				// Process/row/column
		descb[5] = subsize_lappack;  // Local Size
		descb[6] = 0;				// Reserved

	double* wk = new double[lwork];

	//opening the text file to output the results
	ofstream f_out2;
	if (rank == (size/2) -1){
		f_out2.open("Out.txt");
		if (!f_out2.good()) {
				cout << "Error: unable to open output file: output.txt" << endl;
		}
	}

	while (t <= T){
		// Adding the steady state condition
		if (t < T1){
			Fy = -1000.0/T1 * t;
			qy = (-1000.0/T1 * t)*l;
		} else {
			Fy = -1000.0; 
			qy = -1000.0*l;
		}
		
		//Update the force matrix 
		makeForceMat(DoF, F, qy, Fy, subdomain_size, subdomain_start);

		// Copying the original K into Keff as the solver (dpbsv) changes the matrix every iteration
		F77NAME(dcopy)(subdomain_size*6, K, 1, Keff, 1);
		
		//Solving the right hand side of equation 12
		for (int i = 0; i < subdomain_size; ++i){
			F[i] = F[i] + (dt*dt*M[i]*k4[i]);
		}
		
		// Synch all of the procesess before the solver to avoid possible issues
		MPI_Barrier( MPI_COMM_WORLD );
		// Solve the equation using ScaLapack
		F77NAME(pdpbsv)(&U, &DoF, &k, &NRHS, Keff, &ja, desca, F, &ib, descb, wk, &lwork, &INFO);
		if (INFO){cout <<"*****" << INFO << "*****" << endl;}

		// Solve the acceleration and velocity terms note F was overwritten with the values of the displacements
		for (int i = 0; i < subdomain_size; ++i){
			nddot1[i] = k1 * (F[i] - n[i]) - k2 * ndot[i] - (1.0/(2*beta) - 1) * nddot[i];
			ndot[i] = ndot[i] + dt*(1.0 - lambda)*nddot[i] + nddot1[i]*dt*lambda;
		}

		// Copying the results from their corresponding n+1 terms to their n for the next iteration
		F77NAME(dcopy)(subdomain_size, nddot1, 1, nddot, 1);
		F77NAME(dcopy)(subdomain_size, F, 1, n, 1);

		// Solving the parenthesis that mulitplies the mass matrix
		for (int i = 0; i < subdomain_size; ++i){
			k4[i] = k1*n[i] + k2*ndot[i] + k3*nddot[i];
		}

		if (rank == (size/2) - 1){
			// Output to a text file (Only output the term in the middle of the beam)
			f_out2.precision(12);
			f_out2.width(14);
			f_out2 << fixed << n[(DoF)/2-subdomain_start] << " ";
			f_out2 << fixed << t;
			f_out2 << "\n ";
		}
		t = t +dt;
	}
	// Close the file
	if (rank == (size/2) - 1){f_out2.close();}
	Cblacs_gridexit(ctx);
	
	delete[] wk;
}
