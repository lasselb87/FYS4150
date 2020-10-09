//=================================
// Included dependencies
#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <iomanip>
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
//----------------------------------------------------------------------------
{
		if (argc < 5)
		{
				cout << "Bad usage" << endl;
				cout << "Please supply:" << endl;
				cout << "	number of mesh points N"<< endl;
				cout << "	starting point p0"<< endl;
				cout << "	end point pN"<< endl;
				cout << "	oscillator frequencys w"<< endl;
				return 1;
		}
		int N = atoi(argv[1]);
		double p0 = atof(argv[2]);
		double pN = atof(argv[3]);
		double w = atof(argv[4]);

		double h = (pN - p0)/(N+1);
		double d = 2/(h*h);
		double a = -1/(h*h);
		vec x = linspace(p0+h, pN-h, N);  // Interval

		// Declare eigval for storing eigenvalues
		//(values not used, but declaration neccessary)
		vec eigval;
		// Declare eigvec1 for storing non-interacting eigenvectors
		mat eigvec1;
		// Declare eigvec2 for storing interacting eigenvectors
		mat eigvec2;

		// Non-interacting case
		mat A = makeMatrix(d, a, N);
		vec V1 = harmOsc(x, w);
		A.diag(0) += V1;
		jacobiMethod(A, eigval, eigvec1, N);

		// Interacting case
		mat B = makeMatrix(d, a, N);
		vec V2 = coulomb(x, w);
		B.diag(0) += V2;
		jacobiMethod(B, eigval, eigvec2, N);


		// Eigenvectors normalized by dividing them by sqrt(h)
		ofstream myfile;
		myfile.open("./Results/waveFunctions.txt");
		for(int i=0; i<N; i++)
		{
				myfile << x(i) <<
				        " " << eigvec1(i,0)/sqrt(h)<<
				        " " << eigvec2(i,0)/sqrt(h)<< endl;
		}
		myfile.close();

		string command = string("python3 plot.py waveFunctions ") + string(argv[4]);
		system(command.c_str());

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
