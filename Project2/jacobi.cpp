#include "jacobi.h"


//=============================================================================
double maxoffdiag(mat &A, int &k, int &l, int N)
//----------------------------------------------------------------------------
// Find the largest off-diagonal element above the diagonal of the
// given matrix
//
// A - matrix to search
// k,l - indices of largest element
// N - dimension of matrix
//----------------------------------------------------------------------------
{
		double max = 0.0;
		for (int i=0; i<N; i++) {
				for (int j=i+1; j<N; j++) {
						if (fabs(A(i,j)) > max) {
								max = fabs(A(i,j));
								l = i;
								k = j;
						}
				}
		}
		return max;
}

//=============================================================================
void jacobiRotate(mat &A, mat &R, int k, int l, int N)
//----------------------------------------------------------------------------
// Does a single rotation of the matrix A to remove the element A(k,l)
// Also rotates the accompanying set of verctors R
//
// A - matrix to be diagonalized
// R - accompanying set of eigenvectors to be rotated
// k,l - element to be set to zero
// N - dimension of matrix
//----------------------------------------------------------------------------
{
		double a_kk = A(k, k);
		double a_ll = A(l, l);
		double a_kl = A(k, l);
		double t, s, c;
		double tau = (a_ll - a_kk)/(2*a_kl);

		if (tau > 0) {
				t = 1.0/(tau + sqrt(1.0 + tau*tau));
		}
		else {
				t = -1.0/( -tau + sqrt(1.0 + tau*tau));
		}

		c = 1/sqrt(1+t*t);
		s = t*c;

		for (int i=0; i<N; i++) {
				double a_ik_prime = A(i, k);
				double a_il_prime = A(i, l);
				A(k, i) = A(i, k) = c*a_ik_prime - s*a_il_prime;
				A(l, i) = A(i, l) = c*a_il_prime + s*a_ik_prime;

				double r_ik_prime = R(i, k);
				double r_il_prime = R(i, l);
				R(i,k) = c*r_ik_prime - s*r_il_prime;
				R(i,l) = c*r_il_prime + s*r_ik_prime;
		}

		A(k, k) = c*c*a_kk - 2.0*c*s*a_kl + s*s*a_ll;
		A(l, l) = s*s*a_kk + 2.0*c*s*a_kl + c*c*a_ll;
		A(k, l) = A(l, k) = 0.0;

		return;
}

//============================================================================
int jacobiMethod(mat A, vec &eigval, mat &eigvec, int N)
//----------------------------------------------------------------------------
// Diagonalize A and generate a sorted set of eigenvalues and eigenvectors
//
// Input:
// A - matrix to be diagonalized
// eigval - address to put eigenvalues
// eigvec - address to put eigenvectors
// N - dimension of matrix
//
// Return:
// number of similarity transformations
//----------------------------------------------------------------------------
{
		int k, l;
		double max = 1;
		double eps = 1e-10;

		mat R(N,N,fill::zeros);       // For storing unsorted eigenvectors
		R.diag(0).fill(1.);

		int iterations = 0;           // Count number of rotations
		while (max>eps)
		{
				max = maxoffdiag(A, k, l, N);
				jacobiRotate(A, R, k, l, N);
				iterations++;
		}

		vec eigval_R = A.diag(0);   // Unsorted eigenvalues
		eigval = sort(eigval_R);    // Sort eigenvalues
		eigvec.zeros(N, N);

		for(int i = 0; i<N; i++)
		{
				for(int j = 0; j<N; j++)
				{
						if (eigval_R(j) == eigval(i))
						{
								eigvec.col(i) = R.col(j);  // Sort eigenvectors
						}
				}
		}

		return iterations;
}
