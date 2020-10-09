#include <iostream>
#include <armadillo>
#include <chrono>
#include <fstream>

#include "jacobi.h"

using namespace arma;
using namespace std;
using namespace chrono;

int main()
{
		ofstream myfile;
		myfile.open("./Results/benchmark.txt");

		vec eigval;
		mat eigvec;
		int iter;
		int m = 10;
		double average_time_jacobi = 0;
		double average_time_arma = 0;
		// Run Jacobi's method and Armadillo's eig_sym N = 10, 20, ..., 120
		for(int N = 10; N <= 200; N += 10)
		{
				// Run the methods m = 10 times for each N, then average the times
				for(int i = 0; i < m; i++)
				{
						mat A = zeros(N, N);
						for (int j = 0; j < N-1; j++)
						{
								A(j,j) = 2;
								A(j,j+1) = -1;
								A(j+1,j) = -1;
						}
						A(N-1,N-1) = 2;

						auto start = high_resolution_clock::now();
						iter = jacobiMethod(A, eigval, eigvec, N);
						auto finish = high_resolution_clock::now();
						average_time_jacobi += duration<double>(finish - start).count();

						start = high_resolution_clock::now();
						eig_sym(eigval, eigvec, A);
						finish = high_resolution_clock::now();
						average_time_arma += duration<double>(finish - start).count();
				}
				// Average the times
				average_time_jacobi = average_time_jacobi/m;
				average_time_arma = average_time_arma/m;
				// Write to file
				myfile << N << " " << average_time_jacobi << " " << average_time_arma
				       <<" " << iter << endl;
		}
		myfile.close();
		return 0;
}
