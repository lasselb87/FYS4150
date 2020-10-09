//=================================
// Include guard
#ifndef __JACOBI_H_INCLUDED__
#define __JACOBI_H_INCLUDED__

//=================================
// Included dependencies
#include <iostream>
#include <cmath>
#include <armadillo>

using namespace  std;
using namespace  arma;

//=================================
// Function prototypes
double maxoffdiag(mat &A, int &k, int &l, int N);
void jacobiRotate(mat &A, mat &R, int k, int l, int N);
int jacobiMethod(mat A, vec &eigval, mat &eigvec, int N);

#endif // __JACOBI_H_INCLUDED__
