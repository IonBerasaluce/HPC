/* Function Files for task 2
 *	
 * Written by: Ion Berasaluce
 * CID: 00950097
 */

#include "headers.h"
#include <iostream>
#include <fstream>
using namespace std;

#define F77NAME(x) x##_
extern "C" {

	void F77NAME(daxpy)(const int& N, const double& DA, const double* DX, const int& INCX,
		double* DY, const int& INCY);

	void F77NAME(dsbmv)(const char& UPLO, const int& N, const int& K, const double& ALPHA, double* A,
		const int& LDA, double* X, const int& INCX, const double& BETA, double* Y, const int& INCY);

	void F77NAME(dcopy)(const int& N, double* DX, const int& INCX, double* DY,
		const int& INCY);
}

//Explicit Solver for task 2
void explicitSolv(double dt, const int DoF, double* M, double* K,
				  double T, double T1, double l){
/* dt is the timestep selected		
   M is the original mass matrix of the problem
   DoF is the full degrees of freedom of the problem excluding the boundary conditions
   K is the full stiffness matrix of the beam element
   T is the maximum time the load is applied
   T1 is the time at which the steady state begins
   l is the elemental length
*/
	double t = 0; double n[DoF] = {0}; double n_1[DoF] = {0};
	double Fy, qy;
	double F[DoF];
	double sol[DoF];

	// Solve the K-2M parenthesis of the solver equation 
	for (int i = 0; i < DoF; ++i){
		K[6*i + 5] = K[6*i + 5] -2 * M[i];
	}
	

    //opening the text file to output the results
	ofstream f_out2("Out.txt");
	
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
		makeForceMat(DoF, F, qy, Fy, DoF, 0);
      	//solve the right hand side of the equation (9)
      	//n is the vector of displacements at time t, n_1 is the vector of displacements at time t-1
      	//F starts as the Force vector but gets changed into the solution of the rhs of equation 9 
       	F77NAME(dsbmv)('U',DoF, 5, -1.0, K, 6, n, 1, 1.0, F, 1); 
      	F77NAME(dsbmv)('U',DoF, 0, -1.0, M, 1, n_1, 1, 1.0, F, 1); 

   		// divide by the M matrix or multiply but the Inverse mass matrix.
  		for (int i = 0; i < DoF; ++i){
 			sol[i] = F[i] * (1.0/M[i]);
  		}

  		//copying the displacements at t = t to those at t = t-dt and t+dt to t.
		F77NAME(dcopy)(DoF, n, 1, n_1, 1);
  		F77NAME(dcopy)(DoF, sol, 1, n, 1);

      	// Output to a text file (Only output the term in the middle of the beam)
	  	f_out2.precision(12);
	  	f_out2.width(14);
	 	f_out2 << std::fixed << sol[(DoF)/2] << " ";
	  	f_out2 << std::fixed << t;
		f_out2 << "\n ";
		
		//Increase the time counter
		t += dt; 
	} 
	
	//Let the user now the task has finished
	cout << "Written file: Out.txt" << endl;
	// Close file
	f_out2.close();    
}

// Overloaded function for task 2f constains the new arguments min and max
void explicitSolv(double dt, const int DoF, double* M, double* K,
				  double T, double T1, double l, double& max, double& min){
/* dt is the timestep selected		
   M is the original mass matrix of the problem
   DoF is the full degrees of freedom of the problem excluding the boundary conditions
   K is the full stiffness matrix of the beam element
   T is the maximum time the load is applied
   T1 is the time at which the steady state begins
   l is the elemental length
   max and min are the variables containing the maximum and minimum values of the amplitudes
*/
	double t = 0; double n[DoF] = {0}; double n_1[DoF] = {0};
	double Fy, qy;
	double F[DoF];
	double sol[DoF];
	int check = 1;

	// Solve the K-2M parenthesis of the solver equation 
	for (int i = 0; i < DoF; ++i){
		K[6*i + 5] = K[6*i + 5] -2 * M[i];
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
		makeForceMat(DoF, F, qy, Fy, DoF, 0);
      	//solve the right hand side of the equation (9)
      	//n is the vector of displacements at time t, n_1 is the vector of displacements at time t-1
      	//F starts as the Force vector but gets changed into the solution of the rhs of equation 9 
       	F77NAME(dsbmv)('U',DoF, 5, -1.0, K, 6, n, 1, 1.0, F, 1); 
       	// Change for a for-loop
      	F77NAME(dsbmv)('U',DoF, 0, -1.0, M, 1, n_1, 1, 1.0, F, 1); 
    
   		// divide by the M matrix or multiply but the Inverse mass matrix.
  		for (int i = 0; i < DoF; ++i){
 			sol[i] = F[i] * (1.0/M[i]);
  		}

  		//copying the displacements at t = t to those at t = t-dt and t+dt to t.
		F77NAME(dcopy)(DoF, n, 1, n_1, 1);
  		F77NAME(dcopy)(DoF, sol, 1, n, 1);

   		if (t >= T1){
   			// check if we are in the steady state and start to find the minimum and maximum of the oscilations
   			if (check == 1){
  				max = sol[DoF/2];
  				min = sol[DoF/2];
  				check += 1;
  			}
			if (sol[DoF/2] > max){
				max = sol[DoF/2];
			}
			if (sol[DoF/2] < min){
				min = sol[DoF/2];
			}
		}

		//Increase the time counter
		t += dt; 
	} 
}