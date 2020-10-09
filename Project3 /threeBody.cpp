//=================================
// Included dependencies
#include <iostream>
#include <iomanip>
#include <armadillo>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include "Classes/solver.h"
#include "Classes/planet.h"
using namespace std;
using namespace arma;

//=================================
// Global variable(s)
double GM = 4*M_PI*M_PI;      // Gravitational constant*1 solar mass


//===================================
//------------ FUNCTIONS ------------
//===================================

//============================================================================
inline vec newton(vec pos, vec vel)
//----------------------------------------------------------------------------
// Calculate newtonian gravity
//----------------------------------------------------------------------------
{
	double rCube = pow(norm(pos), 3);
	return -GM/rCube*pos;
}
//============================================================================


//============================================================================
//-------------------------------- MAIN --------------------------------------
//============================================================================
int main()
//----------------------------------------------------------------------------
// Simulate the motion of the Earth-Jupiter-Sun system
//----------------------------------------------------------------------------
{

	double T = 22;            // simulation time [years]
	int N = 1000000;          // N - number of integration points
	int sampleN = 10;         // n - sample point interval

	// Initialize celestial bodies in the Solar System
	Planet Earth(vec({1,0,0}), vec({0,2*M_PI,0}), 3.039E-6);
	Planet Jupiter(vec({5,0,0}), vec({0,3,0}), 954.7E-6);
	vector<Planet> solarsystem = vector<Planet>{Earth, Jupiter};
	vector<double> mass{954.7E-6, 10*954.7E-6, 1000*954.7E-6};


	// Solve motion with Verlet's method
	for(int i = 1; i<=3; i++)
	{
		solarsystem[1].M = mass[i-1];
		Solver solver(solarsystem, GM, true);
		solver.solve(2, newton, T, N, sampleN, "./Raw_Data/data" + to_string(i) + ".txt");
	}
	system("python3 plot.py earthAndJupiter");

	// Calculate energy fluctation with Verlet's method
	ofstream myfile;
	myfile.open("./Raw_Data/fluctuation.txt");
	for(int n = 100; n<=1e8; n*=10)
	{
		Solver solver(solarsystem, GM, true);
		solver.solve(2, newton, T, n, n/100, "./Raw_Data/data.txt");
		myfile << n << " " << solver.totalEnergyFluctuation() << "\n";
		cout << n << endl;
	}
	myfile.close();
	system("python3 plot.py fluctuation EarthAndJupiter");

	return 0;
}
//============================================================================
