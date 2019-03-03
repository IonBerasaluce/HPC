/* Function Files for task 4
 *	
 * Written by: Ion Berasaluce
 * CID: 00950097
 */


#include "headers.h"
#include <iostream>
#include <fstream>
#include <mpi.h>
using namespace std;


#define F77NAME(x) x##_
extern "C" {
	void F77NAME(dsbmv)(const char& UPLO, const int& N, const int& K, const double& ALPHA, double* A,
    const int& LDA, double* X, const int& INCX, const double& BETA, double* Y, const int& INCY);

    void F77NAME(dcopy)(const int& N, double* DX, const int& INCX, double* DY,
    const int& INCY);
}

void explicitSolvPar(const int subdomain_size, const int subdomain_start, double* K, double* M, double dt, double T, double T1, int DoF, double l){
/* subdomain size is the domain of the partitioned DoF of the specific process
   subdomain starts is the starting point of the Degrees of Freedom of the process
   DoF is the total number of Degrees of Freedom of the system excluding those at the boundary conditions
   M is the mass matrix of the system
   K is the assembled full stiffness matrix
   dt is the timestep used 
   T is the maximum time 
   T1 is the maximum loading time 
   l is the length of the elements
*/
 	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// Variables required for the solver ie the displacement u at n and the displacement u at n-1
    double t = 0; 						double n[subdomain_size] = {0}; 	
 	double n_1[subdomain_size] = {0};	double F[subdomain_size];
    double sol[subdomain_size];			double transf[3];
    double Fy, qy;					
    
    ofstream f_out2;
    // open teh text file
    if (rank == 0){
      f_out2.open("Out.txt");
      if (!f_out2.good()) {
        cout << "Error: unable to open output file: output.txt" << endl;
      }
    }

    while (t <= T){
	    //adding the steady state condition
	    if (t < T1){
	    	Fy = -1000.0/T1 * t;
	        qy = (-1000.0/T1 * t)*l;
	    } else {
	        Fy = -1000.0; 
	        qy = -1000.0*l;
	    }
	      
	    //make the updated Force vector with the new values of Fy and qy
	    //cout << "Rank: " << rank << "subdomain_size: " << subdomain_size << "subdomain_start: " << subdomain_start << endl;
	    makeForceMat(DoF, F, qy, Fy, subdomain_size, subdomain_start);

	    F77NAME(dsbmv)('U',subdomain_size, 5, -1.0, K, 6, n, 1, 1.0, F, 1);
	    F77NAME(dsbmv)('U',subdomain_size, 0, -1.0, M, 1, n_1, 1, 1.0, F, 1); 

	    // divide by the M matrix or multiply but the Inverse mass matrix.
	    for (int i = 0; i < subdomain_size; ++i){
	          sol[i] = F[i] * (1.0/M[i]);
	    }
	    // Mpi messaging to send the negihbouring elements to the central node from P1 to P0 and vice versa
	    if (rank == 1){
	        for(int i = 0; i < 3; ++i){
	        	// We grab elements 6,7,8 as 0,1,2 are the ones we are getting from rank 0 and 3,4,5 are the central node.
	            transf[i] = sol[i+6];
	        }
	        MPI_Send(&transf, 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	    }

	    if (rank == 0){
	        MPI_Recv(&transf, 3, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        for (int i = 0; i < 3; ++i){
	        	// The values recieved from process 1 are placed in the last 3 locations of process 2
	            sol[subdomain_size - (3 - i)] = transf[i];
	        }  
	    
	        for (int i = 0; i < 3; ++i){
	        	// Again the values required to transfer are the 3rd to last triplet since the last triplet, comes from process 1 and the 
	        	// 2nd to last is the middle node.
	            transf[i] = sol[subdomain_size - (9-i)];
	        }
	        MPI_Send(&transf, 3, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
	    } 

	    if (rank == 1){
	        MPI_Recv(&transf, 3, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        for (int i = 0; i < 3; ++i){
	        	// We place the elements received from process 0 into the first 3 locations of proces 1.
	            sol[i] = transf[i];
	        }   
	    }

		// Copying the displacements at t = t to those at t = t-dt and t+dt to t.
		F77NAME(dcopy)(subdomain_size, n, 1, n_1, 1);
		F77NAME(dcopy)(subdomain_size, sol, 1, n, 1);
	     
	    // Output the results of the displacements of the middle node to a text file
		if (rank == 0){ 
	      f_out2.precision(10);
	      f_out2.width(12);
	      f_out2 << std::fixed << sol[DoF/2] << " ";
	      f_out2 << std::fixed << t << "\n";
	  	}
	      
	    //Increase the time counter
	    t += dt; 
    }
    // Close the file
    if (rank == (size/2) - 1){f_out2.close();} 
 }
