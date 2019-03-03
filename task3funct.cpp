/* Function Files for task 3
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
	void F77NAME(dpbsv)(const char& UPLO, const int& N, const int& KD, const int& NRHS,
		double* AB, const int& LDAB, double* B, const int& ldb, int& INFO);

	void F77NAME(dsbmv)(const char& UPLO, const int& N, const int& K, const double& ALPHA, double* A,
		const int& LDA, double* X, const int& INCX, const double& BETA, double* Y, const int& INCY);

	void F77NAME(dcopy)(const int& N, double* DX, const int& INCX, double* DY,
		const int& INCY);
}

// Implicit solver for task 3
void implicitSolv(double dt, const int DoF, double* M, double* K, double T, double T1, double l){
/* dt is the solver timestep 
   DoF is the total number of DoF of the system excluding those at the boundary conditions
   M is the mass matrix of the system
   K is the assembled full stiffness matrix
   T is the maximum time the load is applied
   T1 is the time at which the steady state begins
   l is the elemental length
*/
	// defining the specific variables required for the solver.
	double beta = 0.25; double lambda = 0.5; double t = 0; int INFO; double Keff[6*(DoF)];
	
	double n[DoF] = {0};double ndot[DoF] = {0}; 
	double nddot1[DoF] = {0};  double nddot[DoF] = {0};
	
	double k1 = 1.0/(beta*dt*dt); double k2 = 1.0/(dt*beta); double k3 = (1.0/(2*beta) - 1); double k4[DoF] = {0}; double F[DoF];
	double Fy,qy;


	//Make Keff into K
	for (int i = 0; i < DoF; ++i){
		for (int j = 0; j < 6; ++j){
			if (j == 5){ 
				K[6*i + j] =  K[6*i + j] + (1.0/beta) * M[i];
			}
		}
	}

	//opening the text file to output the results
	ofstream f_out3("Out.txt");
	
	while (t <= T){
		
		if (t < T1){
			Fy = -1000.0/T1 * t;
			qy = (-1000.0/T1 * t)*l;
		} else {
			Fy = -1000.0; 
			qy = -1000.0*l;
		}
		
		//Update the force atrix 
		makeForceMat(DoF, F, qy, Fy, DoF, 0);

		// Copying the original K into Keff as the solver (dpbsv) changes the matrix every iteration
		F77NAME(dcopy)((DoF)*6, K, 1, Keff, 1);

		//Solving the right hand side of equation 12
		for (int i = 0; i < DoF; ++i){
			F[i] = F[i] + (dt*dt*M[i]*k4[i]);
		}

		// Solve Equation 12
		F77NAME(dpbsv)('U', DoF, 5, 1, Keff, 6, F, DoF, INFO);
  		if (INFO){cout <<"*****" << INFO << "*****" << endl;}
  		// Solve the acceleration and velocity terms at that timestep
		for (int i = 0; i < DoF; ++i){
			nddot1[i] = k1 * (F[i] - n[i]) - k2 * ndot[i] - (1.0/(2*beta) - 1) * nddot[i];
			ndot[i] = ndot[i] + dt*(1.0 - lambda)*nddot[i] + nddot1[i]*dt*lambda;
		}

		//copying the results from their corresponding n+1 terms to their n for the next iteration
		F77NAME(dcopy)(DoF, nddot1, 1, nddot, 1);
		F77NAME(dcopy)(DoF, F, 1, n, 1);

		// solving the parenthesis that multiplies the mass matrix
		for (int i = 0; i < DoF; ++i){
			k4[i] = k1*n[i] + k2*ndot[i] + k3*nddot[i];
		}
		
		// Output to a text file (Only output the term in the middle of the beam)
	  	f_out3.precision(12);
	  	f_out3.width(14);
	 	f_out3 << fixed << n[(DoF)/2] << " ";
	  	f_out3 << fixed << t;
		f_out3 << "\n "; 
		
		// Increase the time counter
		t += dt; 
	}
	// Let the user now the task has finished
	cout << "Written file: Out.txt" << endl;
	// Close file
	f_out3.close();
}

// Overloaded function for task 3c
void implicitSolv(double dt, const int DoF, double* M, double* K, double T, double T1, double l, 
				  double& max, double& min){
/* dt is the solver timestep 
   DoF is the total number of DoF of the system including those at the boundary conditions
   M is the mass matrix of the system
   K is the assembled full stiffness matrix
   T is the maximum time the load is applied
   T1 is the time at which the steady state begins
   l is the elemental length
   max and min to compute the amplitudes of the sytem
*/
	// defining the specific variables required for the solver.
	double beta = 0.25; double lambda = 0.5; double t = 0; int INFO; double Keff[6*(DoF)];
	
	double n[DoF] = {0};double ndot[DoF] = {0}; 
	double nddot1[DoF] = {0};  double nddot[DoF] = {0};
	
	double k1 = 1.0/(beta*dt*dt); double k2 = 1.0/(dt*beta); double k3 = (1.0/(2*beta) - 1); double k4[DoF] = {0}; double F[DoF];
	double Fy,qy;
	int check = 1;

	//Make Keff into K
	for (int i = 0; i < DoF; ++i){
		for (int j = 0; j < 6; ++j){
			if (j == 5){ 
				K[6*i + j] =  K[6*i + j] + (1.0/beta) * M[i];
			}
		}
	}

	while (t <= T){
		
		if (t < T1){
			Fy = -1000.0/T1 * t;
			qy = (-1000.0/T1 * t)*l;
		} else {
			Fy = -1000.0; 
			qy = -1000.0*l;
		}
		
		//Update the force atrix 
		makeForceMat(DoF, F, qy, Fy, DoF, 0);

		// Copying the original K into Keff as the solver (dpbsv) changes the matrix every iteration
		F77NAME(dcopy)((DoF)*6, K, 1, Keff, 1);

		//Solving the right hand side of equation 12
		for (int i = 0; i < DoF; ++i){
			F[i] = F[i] + (dt*dt*M[i]*k4[i]);
		}

		// Solve Equation 12 of the handout
		F77NAME(dpbsv)('U', DoF, 5, 1, Keff, 6, F, DoF, INFO);
  		if (INFO){cout <<"*****" << INFO << "*****" << endl;}
  		// Solve the acceleration and velocity terms at that timestep
		for (int i = 0; i < DoF; ++i){
			nddot1[i] = k1 * (F[i] - n[i]) - k2 * ndot[i] - (1.0/(2*beta) - 1) * nddot[i];
			ndot[i] = ndot[i] + dt*(1.0 - lambda)*nddot[i] + nddot1[i]*dt*lambda;
		}

		//copying the results from their corresponding n+1 terms to their n for the next iteration
		F77NAME(dcopy)(DoF, nddot1, 1, nddot, 1);
		F77NAME(dcopy)(DoF, F, 1, n, 1);

		// Solving the parenthesis that multiplies the mass matrix in equation 12 for the next iteration
		for (int i = 0; i < DoF; ++i){
			k4[i] = k1*n[i] + k2*ndot[i] + k3*nddot[i];
		}
		
		if (t >= T1){
			// Just like in 2f we check to see if we are in the steady response of the system
   			if (check == 1){
  				max = n[DoF/2];
  				min = n[DoF/2];
  				check += 1;
  			}
			if (n[DoF/2] > max){
				max = n[DoF/2];
			}
			if (n[DoF/2] < min){
				min = n[DoF/2];
			}
		}

		// Increase the time counter
		t += dt; 
	}
}