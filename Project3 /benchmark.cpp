#include <iostream>
#include <iomanip>
#include <armadillo>
#include <fstream>
#include <cmath>
#include <string>
#include <chrono>
using namespace std;
using namespace arma;

//=================================
// Global variable(s)
double G = 4*M_PI*M_PI;       // Gravitational constant


//=================================
// Function prototypes
vec acceleration(vec pos);
void euler(vec pos, vec vel, vec acc(vec), int N, double dt);
void verlet(vec pos, vec vel, vec acc(vec), int N, double dt);


//============================================================================
//-------------------------------- MAIN --------------------------------------
//============================================================================
int main()
//----------------------------------------------------------------------------
// Benchmark the CPU time of both Euler's and Verlet's method
//----------------------------------------------------------------------------
{

	double T = 5; // Simulation time
	ofstream myfile;
	myfile.open("./Raw_Data/benchmark.txt");

	// Solve motion and benchmark
	for(int N = 10; N <= 10000000; N *= 10)
	{
		double dt = T/N;    // Time step

		// Initial conditions
		vec pos_euler({1, 0, 0});
		vec vel_euler({0, 2*M_PI, 0});
		vec pos_verlet({1, 0, 0});
		vec vel_verlet({0, 2*M_PI, 0});

		// Benchmark Euler
		auto start1 = chrono::high_resolution_clock::now();
		euler(pos_euler, vel_euler, acceleration, N, dt);
		auto finish1 = chrono::high_resolution_clock::now();
		chrono::duration<double> elapsed1 = finish1 - start1;

		// Benchmark Verlet
		auto start2 = chrono::high_resolution_clock::now();
		verlet(pos_verlet, vel_verlet, acceleration, N, dt);
		auto finish2 = chrono::high_resolution_clock::now();
		chrono::duration<double> elapsed2 = finish2 - start2;

		// Write to file
		myfile << N << " " << elapsed1.count() << " " << elapsed2.count() << "\n";
	}
	myfile.close();
	system("python3 plot.py benchmark");
}
//============================================================================


//===================================
//------------ FUNCTIONS ------------
//===================================

//============================================================================
vec acceleration(vec pos)
//----------------------------------------------------------------------------
// Calculate newtonian gravity
//----------------------------------------------------------------------------
{
	double M = 1;
	double rCube = pow(norm(pos), 3);
	return -G*M/rCube*pos;
}
//============================================================================


//============================================================================
void euler(vec pos, vec vel, vec acc(vec), int N, double dt)
//----------------------------------------------------------------------------
// Euler's method
//----------------------------------------------------------------------------
{
	for(int i=0; i<N-1; i++)
	{
		vel = vel + acc(pos)*dt;
		pos = pos + vel*dt;
	}
}
//============================================================================


//============================================================================
void verlet(vec pos, vec vel, vec acc(vec), int N, double dt)
//----------------------------------------------------------------------------
// Verlet's method
//----------------------------------------------------------------------------
{
	vec prevPos(3,fill::zeros);
	for(int i=0; i<N-1; i++)
	{
		prevPos = pos;
		pos = pos + vel*dt + 0.5*acc(pos)*dt*dt;
		vel = vel + 0.5*(acc(pos) + acc(prevPos))*dt;
	}
}
//============================================================================
