/* Header file for the HPC coursework functions 
 *
 *  By: Ion Berasaluce 
 *  CID: 00950097
 *  Due Date: 26 March
 */

//
void makeStiffMat(int DoF, double* K, double A, double E, double I, double l);
//
void makeForceMat(int DoF, double* F,double qy, double Fy, int subdomain_size, int subdomain_start);
//
void makeMassMat(int DoF, double* M, double l, double dt, double rho, double A);
//
void invertMMatrix(double* M, int DoF);
//
void explicitSolv(double dt, const int DoF, double* M, double* K, double T, double T1,
				  double l);
//
void explicitSolv( double dt, const int DoF, double* M, double* K,
				  double T, double T1, double l, double& max, double& min);
//
void implicitSolv(double dt, const int DoF, double* M, double* K, double T, double T1, double l);
//
void implicitSolv(double dt, const int DoF, double* M, double* K, double T, double T1, double l,
				  double& max, double& min);
//
void decompose_domain(const int domain_size, int world_rank, int world_size, int& subdomain_size, 
					  int& subdomain_start);
//
void ExportFile(double dataX, double dataY);
//
void printMatrix(int DoF, double* A);
//
void explicitSolvPar(const int subdomain_size, const int subdomain_start, double* K, double* M, 
					 double dt, double T, double T1, int DoF, double l);
//
void implicitSolvPar(const int subdomain_size, const int subdomain_start, double* K, double* M,
					 double dt, double T, double T1, double beta, const int DoF, double l);
//
void printVector(double* V, int DoF);
//
void zeros(double* M, int size);