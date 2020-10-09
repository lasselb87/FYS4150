//=================================
// Included dependencies
#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "jacobi.h"

using namespace  std;
using namespace  arma;

//=================================
// Function prototypes
mat makeMatrix(double d, double a, int N);
vec harmOsc(vec x, double w);
vec coulomb(vec x, double w);


// Begin main program
//=============================================================================
int main(int argc, char *argv[])
//----------------------------------------------------------------------------
// Solve the quantum dot eigenvalue problem using Jacobi's method
//
// Command line arguments:
// N - number of mesh points
// p0 - interval starting point
// pN - interval end point
// w - oscillator frequency
// Which potential: 1 for harmOsc, 2 for coulomb
//----------------------------------------------------------------------------
{
	//Checks if the user is supplying enough commandline arguments
	if (argc < 6)
	{
		cout << "Bad usage" << endl;
		cout << "Please supply:" << endl;
		cout << "Number of mesh points N"<< endl;
		cout << "Interval starting point p0"<< endl;
		cout << "Interval	end point pN"<< endl;
		cout << "Oscillator frequency w"<< endl;
		cout << "Which potential: 1 for harmOsc, 2 for coulomb"<< endl;
		return 1;
	}
	int N = atoi(argv[1]);
	double p0 = atof(argv[2]);
	double pN = atof(argv[3]);
	double w = atof(argv[4]);
	int typePotential = atoi(argv[5]);

	double h = (pN - p0)/(N+1);
	double d = 2/(h*h);
	double a = -1/(h*h);
	vec x = linspace(p0+h, pN-h, N);          // Interval

	vec potential;
	if (typePotential == 1)
	{
		potential = harmOsc(x,w);
	}
	else if (typePotential == 2)
	{
		potential = coulomb(x,w);
	}
	else
	{
		cout << "Bad usage" << endl;
		cout << "Please use 1 for harmOsc, 2 for coulomb" << endl;
		return 1;
	}

	vec eigval;                   // Declare eigval for storing eigenvalues
	mat eigvec;                   // Declare eigvec for storing eigenvectors

	mat A = makeMatrix(d, a, N);
	A.diag(0) += potential;

	// Call the method to yield eigenvalues and eigenvectors
	int iter = jacobiMethod(A, eigval, eigvec, N);

	cout << "Number of mesh points: " << N << endl;
	cout << "Number of similarity transformations: " << iter << endl;
	cout << setprecision(7) << "Ground state eigval: " << eigval(0) << endl;
	cout << setprecision(7) << "1st excited state eigval: " << eigval(1) << endl;
	cout << setprecision(7) << "2nd excited state eigval: " << eigval(2) << endl;
	cout << setprecision(7) << "3rd excited state eigval: " << eigval(3) << endl;
	return 0;
}
// End main program

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


//=============================================================================
vec harmOsc(vec x, double w)
// Harmonic oscillator potential describing one electron
// (or twoNon-interacting electons)
//
// x - interval the potential is defined
// w - oscillator frequency
//----------------------------------------------------------------------------
{
	return w*w*x%x;
}
// End function

//=============================================================================
vec coulomb(vec x, double w)
// Describes twi interaction electrons in a harmonic
//
// x - interval the potential is defined
// w - oscillator frequency
//----------------------------------------------------------------------------
{
	return w*w*x%x + 1/x;
}
// End function
