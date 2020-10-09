//=================================
// Included dependencies
#include <iostream>
#include <cmath>
#include <armadillo>
#include <iomanip>
#include "jacobi.h"

using namespace  std;
using namespace  arma;


int main()
{
		// Test maxoffdiag
		bool passed = true;

		int n = 4;
		mat A = zeros(n, n);
		for (int i = 0; i < n; i++)
		{
				for (int j = 0; j < n; j++)
				{
						A(i,j) = 10*i + j;
				}
		}

		int k = 0;
		int l = 0;

		double max = maxoffdiag(A, k, l, n);
		if (k != 3 or l != 2) {passed = false;}

		A(1, 2) = -100;
		max = 0;

		max = maxoffdiag(A, k, l, n);
		if (k != 2 or l != 1) {passed = false;}

		if (passed) {cout << "Passed: Found larget element" << endl;}
		else{cout << "Failed: Didn't find larget element" << endl;}

		// Test eigenvalues
		passed = true;

		n = 4;
		A = zeros(n, n);
		for (int i = 0; i < n-1; i++)
		{
				A(i,i) = 2;
				A(i,i+1) = -1;
				A(i+1,i) = -1;
		}
		A(n-1,n-1) = 2;

		mat eigvec;
		vec eigval;
		jacobiMethod(A, eigval, eigvec, n);

		double eps = 1e-6;
		double analytical;
		for (int i = 0; i < n; i++)
		{
				analytical = 2*(1 - cos((i+1)*M_PI/(n+1)));
				if (abs(eigval(i) - analytical)> eps)
				{
						passed = false;
				}
				cout << "Numerical: "  << setprecision(4) << eigval(i)
				     << "  Analytical: " << setprecision(4) << analytical
				     << "  Rel error: " << setprecision(4)
				     << abs(analytical - eigval(i))/(analytical) << endl;
		}
		if (passed) {cout << "Passed: All eigenvalues are correct" << endl;}
		else{cout << "Failed: Incorrect eigenvalues" << endl;}

		// Test orthonormality
		passed = true;

		A =  eigvec.t()*eigvec;
		for (int i = 0; i < n; i++)
		{
				for (int j = 0; j < n; j++)
				{
						if (i==j)
						{
								if (abs(A(i,j) - 1)>eps) {passed = false;}
						}
						else
						{
								if (abs(A(i,j))>eps) {passed = false;}
						}
				}
		}
		if (passed == true)
		{
				cout << "Passed: Eigenvectors are orthonormal" << endl;
		}
		else
		{
				cout << "Failed: Eigenvectors are not orthonormal" << endl;
		}

		return 0;
}
