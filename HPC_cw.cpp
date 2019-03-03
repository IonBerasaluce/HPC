/* HPC Finite Element Beam Coursework 
 *
 *  By: Ion Berasaluce 
 *  CID: 00950097
 *  Due Date: 26 March
 */

#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <fstream>
#include <array>
using namespace std;

//include the function heasders required to run the code
#include "headers.h"
#include <mpi.h>

#include <boost/program_options.hpp>

// An alias to reduce typing
namespace po = boost::program_options;

// Declaring the LAPACK routines used
#define F77NAME(x) x##_
extern "C" {
	void F77NAME(dpbsv)(const char& UPLO, const int& N, const int& KD, const int& NRHS,
		double* AB, const int& LDAB, double* B, const int& ldb, int& INFO);

  void F77NAME(dsbmv)(const char& UPLO, const int& N, const int& K, const double& ALPHA, double* A,
    const int& LDA, double* X, const int& INCX, const double& BETA, double* Y, const int& INCY);

  void F77NAME(dcopy)(const int& N, double* DX, const int& INCX, double* DY,
    const int& INCY);

  void F77NAME(pdpbsv)(const char& UPLO, const int& N, const int& nrhs, double* a, const int&  ja, const int* desc_a,
                        double* b, const int& ib, const int* desc_b, double* work, const int& lwork, int& INFO);

}

// main
int main(int argc, char *argv[]){
	po::options_description desc("Finite Element Beam Equation Solver");
	desc.add_options()
	("length",  po::value<double>(),    "Length of the cantilever.")
	("no_elements",  po::value<int>(),    "Number of elements to discretize.")
	("Area",  po::value<double>(),    "Cross-section Area.")
	("MomInertia",  po::value<double>(),    "Moment of Inertia of the cross-section.")
	("Emodulus",  po::value<double>(),    "Youngs Modulus.")
	("density",  po::value<double>(),    "Density.")
	("Max_Time",  po::value<double>(),    "Max-Time.")
	("Interval",  po::value<int>(),    "Interval.")
	("T1",  po::value<double>(),    "Time at which steady state starts.")
	("Static",                     "Solve Static or Dynamic Equations.")
  ("Explicit",                     "Solve the Dynamic Equations using an Explicit Solver.")
  ("Implicit",                     "Solve the Dynamic Equations using an Implicit Solver.")
  ("Serial",                     "Solve the task in serial.")
  ("Parallel",                     "Solve the task in parallel.")
  ("Task2-f",                     "Compute the amplitude of the oscillations for the Explicit Solver.")
	("Task3-c",                     "Compute the amplitude of the oscillations for the Implicit Solver.")
	("help",                       "Print help message.");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

// Help output
	if (vm.count("help")){
		cout << desc << endl;
		return 0;
	}

// Assign the labels to the variables and
// if no numbers are inputted use the default nubers
	const double L      = vm.count("length")      ? vm["length"]       .as<double>()  : 10;
	      int    N 	    = vm.count("no_elements") ? vm["no_elements"]  .as<int>()     : 24;
	const double A 	    = vm.count("Area")        ? vm["Area"]         .as<double>()  : 0.012;
	const double I 	    = vm.count("MomInertia")  ? vm["MomInertia"]   .as<double>()  : 0.0000144;
  const double E      = vm.count("Emodulus")    ? vm["Emodulus"]     .as<double>()  : 210000000000; //GPa
  const double rho    = vm.count("density")     ? vm["density"]      .as<double>()  : 7850; //kg/m^3
  const double T      = vm.count("Max_Time")    ? vm["Max_Time"]     .as<double>()  : 3; //time in secs
  const double T1     = vm.count("T1")    		  ? vm["T1"]     		   .as<double>()  : 1; //time in sece
        int    Nt 	  = vm.count("Interval")    ? vm["Interval"]     .as<int>()     : 100000; //

// negatives check
  if (L < 0 || N < 0 || A < 0 || I < 0 || E < 0 || rho < 0 || T < 0 || Nt < 0){
  	cout << "Cannot accept any negative inputs" << endl;
  	return 1;
  }

// load over a node check (also ensures more than 1 element in use for the parallel programming)
  if (N%2 != 0){
  	cout << "The number of elements must be even to ensure the point load is over a node" << endl;
  	return 1;
  }

// Obtain the distributed load, Point Load
  double Fy = -1000;
  double qy = -1000;

// Calculating the Degrees of Freedom (without boundary conditions)
  const int DoF = ((N + 1) * 3)-6; 
  double l = L/(float)N;

  double K[DoF*6]; 
  double F[DoF];

// Make the Force and Stiffness Matrices 
  //
  makeStiffMat(DoF, K, A, E, I, l);

// Solve the static problem (Task1)
  if (vm.count("Static")){
  	cout << "Solving the Static problem" << endl;
  	int INFO;
    // Creating and assembling the elemental and global force vector for the problem
 	  makeForceMat(DoF, F, qy*l, Fy, DoF, 0);
    // Using the Lapack solver to obtain the solution to Ku = F problem (K is triangular banded symmetric matrix)
    // and F is a vector
    F77NAME(dpbsv)('U', DoF, 5, 1, K, 6, F, DoF, INFO);
    if (INFO){cout <<"*****" << INFO << "*****" << endl;}

  // Output to a text file
  	ofstream f_out1("output.txt");

  	if (!f_out1.good()) {
  		cout << "Error: unable to open output file: output.txt" << endl;
  	} else {	
  		double l1 = 0;
  		for (int i = 1; i < DoF; i+=3) {
  			l1 = l1 + l;
  			f_out1.precision(10);
  			f_out1.width(12);
  			f_out1 << fixed << F[i];
  			f_out1.precision(4);
  			f_out1.width(7);
  			f_out1 << fixed <<l1 << " " << endl;
  		}

      // Close file
  		f_out1.close();

  		cout << "Written file: output.txt" << endl;
  	}
  	return 0;
  }
  
  if (T1 > T){
    cout << "Transient response time greater than total time!" << endl;
    return -1;
  }

// Solve the problem using an Explicit solver (Task2)
 	if (vm.count("Explicit") && vm.count("Serial")){
		cout << "Solving using an Explicit solver" << endl;

		// Make Mass matrix twice as 1 is used for inverse(Mcopy) and the other by itself (M)
		double dt = T/Nt;
		double M[DoF];
		makeMassMat(DoF, M, l, dt, rho, A);

		// Solve using Explicit integration scheme
		explicitSolv(dt, DoF, M, K, T, T1, l);

		return 0;
	}

  // Find the amplitudes of the Oscillations at different times (Task2 f)
  if (vm.count("Task2-f")){
    double M[DoF];
    double max, min;
    //opening the text file to output the results
    ofstream f_out("Out.txt");

    for (int i = 1; i < 66; ++i){
      int Ntf = i*Nt;
      double dt = ((double) i)/((double) Ntf);
      makeMassMat(DoF, M, l, dt, rho, A);
      makeStiffMat(DoF, K, A, E, I, l);

      max = 0; min = 0;

      explicitSolv(dt, DoF, M, K, (double) i, ((double)i)/2.0, l, max, min);
      // Output the max and min to a text file, the amplitude is then processed in python
      if (!f_out.good()) {
        cout << "Error: unable to open output file: output.txt" << endl;
      } else {  
        f_out.precision(12);
        f_out.width(14);
        f_out << fixed << max << " ";
        f_out << fixed << min << " ";
        f_out << fixed << i;
        f_out << "\n ";
      }
    }
    //Let the user now the task has finished
    cout << "Written file: Out.txt" << endl;
    // Close file
    f_out.close(); 
    return 0; 
  }

// Solve the problem using an Implicit solver (Task3)
  if (vm.count("Implicit") && vm.count("Serial")){
	  cout << "Solving using an Implicit solver" << endl;
	  
    double dt = T/Nt;
    double M[DoF];
    makeMassMat(DoF, M, l, dt, rho, A);
    // Solve using an explicit solver 
    implicitSolv(dt, DoF, M, K, T, T1, l);
    return 0;
  }

  if (vm.count("Task3-c")){
    ofstream f_out("Out.txt");
    double max, min;
    double M[DoF];
    
    // Iterating through multiple Tl to obtain the variation in amplitude
    for (int i = 1; i < 66; ++i){
      int Ntf = i*Nt;
      double dt = ((double) i)/((double) Ntf);
      makeMassMat(DoF, M, l, dt, rho, A);
      max = 0; min = 0;

      implicitSolv(dt, DoF, M, K, (double)i, ((double)i)/2.0, l, max, min);

      // Output to a text file (Only output the term in the middle of the beam)
      f_out.precision(12);
      f_out.width(14);
      f_out << fixed << max << " ";
      f_out << fixed << min << " ";
      f_out << fixed << i;
      f_out << "\n "; 
    }

    //Let the user now the task has finished
    cout << "Written file: Out.txt" << endl;
    // Close file
    f_out.close();
    return 0;
  }

// Solve the problem using an Explicit solver in parallel (Task4)
  if (vm.count("Explicit") && vm.count("Parallel")){

    // Common variables to both processes.
    int size = 0;
    int rank = 0;
    double dt = T/Nt;

    // Initialize the MPI interface
    int err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS) {
      cout << "Failed to initialise MPI" << endl;
      return -1;
    }

    // Get the rank and comm size on each process.
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // We split the domain into 2 equal subdomains by repeating the middle node.
    int subdomain_start, subdomain_size;
    subdomain_size = N/size + 1;
    subdomain_start = subdomain_size * rank;
    
    // Ensure that if we are in the second process we also contain the middle node.
    if (rank >= size/2){
      subdomain_start = subdomain_start - 3;
    }

    subdomain_size = subdomain_size*3;
    subdomain_start = subdomain_start*3;

    double M[subdomain_size];
    double Mcopy[subdomain_size];
    double K[subdomain_size*6];
    
    // Compute the local mass, force and stiffness matrix
    makeMassMat(subdomain_size, M, l, dt, rho, A);
    makeMassMat(subdomain_size, Mcopy, l, dt, rho, A);
    makeStiffMat(subdomain_size, K, A, E, I, l);

    // Solve the K-2M parenthesis of the solver equation 
    for (int i = 0; i < subdomain_size; ++i){
        K[6*i + 5] = K[6*i + 5] -2 * M[i];
    }
    
    // Implement the parallel explicit solver which makes the output directly to a .txt file
    explicitSolvPar(subdomain_size, subdomain_start, K, M, dt, T, T1, DoF, l);
    
    MPI_Finalize();
    return 0;
  }
    
// Solve the problem using an Implicit solver in parallel(Task5)

  if (vm.count("Implicit") && vm.count("Parallel")){
    // Variable definition for this problem
    int rank, size;
    double dt = T/Nt;
    double M[DoF];
    double beta = 0.25;
    makeStiffMat(DoF, K, A, E, I , l);
    makeMassMat(DoF, M, l, dt, rho, A);

    //Make Keff into K, solver will then need to copy the corresponding terms
    for (int i = 0; i < DoF; ++i){
      for (int j = 0; j < 6; ++j){
        if (j == 5){ 
          K[6*i + j] =  K[6*i + j] + (1.0/beta) * M[i];
        }
      }
    }

    // Initialize the MPI interface 
    int err = MPI_Init(&argc, &argv);
      if (err != MPI_SUCCESS) {
        cout << "Failed to initialise MPI" << endl;
      return -1;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // Create the subdivision of the degrees of freedom between the 2 processes with no repeated nodes in this case
    int subdomain_size = (float)N/(float)size;
    int subdomain_start = subdomain_size * rank;

    if (rank == (size -1)){
      subdomain_size = (N - 1) - subdomain_size*rank;
    }

    subdomain_size = subdomain_size*3;
    subdomain_start = subdomain_start*3;

    // Call the implicit solver function
    implicitSolvPar(subdomain_size, subdomain_start, K, M, dt, T, T1, beta, DoF, l);

    cout << "Complete" << endl;

    MPI_Finalize();
    return 0;
  }

}


    


