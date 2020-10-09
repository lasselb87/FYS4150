//=================================
// Included dependencies
#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <chrono>
#include "jacobi.h"

using namespace  std;
using namespace  arma;

//=================================
// Function prototypes
vec analyticalEigval(double d, double a, int N);
mat makeMatrix(double d, double a, int N);

// Begin main program
//=============================================================================
int main(int argc, char *argv[])
//----------------------------------------------------------------------------
// Solve the buckling beam eigenvalue problem using Jacobi's method
//
// Command line arguments:
// N - number of mesh points
//----------------------------------------------------------------------------
{
		if (argc < 2)
		{
				cout << "Bad usage:" << endl;
				cout << "Supply number of mesh points N on cmd line"<< endl;
				return 1;
		}
		int N = atoi(argv[1]);

		double rho_0 = 0.0;
		double rho_N = 1.0;
		double h = (rho_N-rho_0)/(N+1);
		double d = 2/(h*h);
		double a = -1/(h*h);
		vec eigval;           // Declare eigval for storing eigenvalues
		mat eigvec;           // Declare eigvec for storing eigenvectors

		mat A = makeMatrix(d, a, N);

		auto start = chrono::high_resolution_clock::now();
		int iter = jacobiMethod(A, eigval, eigvec, N);
		auto finish = chrono::high_resolution_clock::now();
		chrono::duration<double> elapsed = finish - start;

		vec ana_eigval = analyticalEigval(d, a, N);

		cout << "Number of mesh points: " << N << endl;
		cout << "Number of similarity transformations: " << iter << endl;
		cout << "Elapsed time (s): " << elapsed.count() << endl;
		cout << "Eigenvalues:" << endl;
		cout << "Numerical" << " " << "Analytical" << endl;
		cout << eigval(0) << " " << ana_eigval(0) << endl;
		cout << eigval(1) << " " << ana_eigval(1) << endl;
		cout << eigval(2) << " " << ana_eigval(2) << endl;

		return 0;
}
// End main program

//=============================================================================
vec analyticalEigval(double d, double a, int N)
//----------------------------------------------------------------------------
// Compute analytical eigenvalues of the buckling beam problem
//
// d - diagonal elements
// a - off-diagonal elements
// N - number of mesh points
//----------------------------------------------------------------------------
{
		vec ana_eigval(N);
		for (int j=1; j<N+1; j++) {
				ana_eigval[j-1] = d+2*a*cos((j*M_PI)/(N+1));
		}
		return ana_eigval;
}
// End function

//=============================================================================
mat makeMatrix(double d, double a, int N)
//----------------------------------------------------------------------------
// Construct a tridiagonal matrix
//
// d - diagonal elements
// a - off-diagonal elements
// N - number of mesh points -> determines matrix dimensionality
//----------------------------------------------------------------------------
{
		mat A(N,N, fill::zeros);
		A.diag(-1).fill(a);
		A.diag(0).fill(d);
		A.diag(1).fill(a);
		return A;
}
// End function
