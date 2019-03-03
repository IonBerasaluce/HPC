//base function file ie. functions that must compile every time the code runs

#include "headers.h"
#include <iostream>
using namespace std;

// Assemble the Stiffnes Matrix with DoF degrees of freedom A area, E youngs Modulus, I moment of Inertia and 
// elemental lenght l.
void makeStiffMat(int DoF, double* K, double A, double E, double I, double l) { 
	// Defining and allocating the different submatrices of the assembled stiff matrix
	double mat2[9] = {-(A * E)/l, 0, 0, 0, -(12 * E * I)/(l*l*l), -(6 * E * I)/(l*l), 0, (6 * E * I)/(l*l), (2 * E * I)/(l)};
	double mat4[6] = {2*(A * E)/l, 0, (24 * E * I)/(l*l*l), 0, 0, (8 * E * I)/(l)};

	int it2 = 9; int it4 = 0;
	for (int i = 0; i < DoF; i++){
		for (int j = 0; j < 6; j++){ 
			if (j >= 5 - (i%3)){
				if (i%3 == 0) {it4 = 0;}
                // We are in the bottom triangle we use mat4
				K[6*i + j] = mat4[it4];
				++it4;
			} else if (j >= 2 - (i%3)){
                // We are in the bottom triangle of the top half of the matrix use mat2
				if (i%3 == 0 && it2 == 9) {it2 = 0;}
				K[6*i + j] = mat2[it2];
				++it2;
			} else {
                // The remaining terms are 0 
				K[6*i + j] = 0;
			}   
		}
	}
}

// Make Mass Matrix into array M of length DoF, element length l, timestep dt, area A and density rho 
void makeMassMat(int DoF, double* M, double l, double dt, double rho, double A){
	double k1 = rho * A * l;
	double alpha = 1.0/24.0;
	int ii = -3;

	for (int i = 0; i < DoF; ++i){
		if (i%3 == 0){
			ii += 3;
		}
		if ((i-ii) < 2){
			M[i] = k1/(dt*dt);
		} else {
			M[i] = 2*k1*alpha*l*l/(dt*dt);
		}
	}
}

// Make force matrix of length subdomain_size, starting at subdomain_start and and add a point load in the central node
// of a problem with DoF degrees of freedom.
void makeForceMat(int DoF, double* F,double qy, double Fy, int subdomain_size, int subdomain_start){
	double mat2[3] = {0, qy, 0};
	int itr1 = 0;
    for (int i = subdomain_start; i < subdomain_start+subdomain_size; i++) {
        if ((i-subdomain_start)%3 == 0) {itr1 = 0;}
	F[i-subdomain_start] = mat2[itr1];
	itr1++;
     // Adding the Point Load to the middle of the beam
	    if (i == (DoF/2)){
	        F[i-subdomain_start] += Fy;
	   	}
    }
}

// Print column-major, banded matrix with DoF columns 
void printMatrix(int DoF, double* A) {
	for (int i = 0; i < 6; i++) {
		cout << "[ ";
		for (int j = 0; j < DoF; j++) {
			cout << A[j * 6 + i] << "  ";
		}
		cout << " ]" << endl;
	}
	cout << "\n" << endl;
}

// Print Vector V of length DoF
void printVector(double* V, int DoF){
	for (int i = 0; i < DoF; ++i)
	{
		cout << V[i] << " ";
	}
	cout << "\n";
}

// Set an array of size size to all zeros.
void zeros(double* M, int size){
	for (int i = 0; i < size; ++i){
		M[i] = 0;
	}
}


