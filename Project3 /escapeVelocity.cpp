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
#include "./Classes/solver.h"
#include "./Classes/planet.h"
using namespace std;
using namespace arma;

//=================================
// Global variable(s)
double G = 4*M_PI*M_PI;      // Gravitational constant


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
	return -G/rCube*pos;
}
//============================================================================

//============================================================================
inline vec newtonBeta1(vec pos, vec vel)
//----------------------------------------------------------------------------
// Calculate inverse-square -> inverse-cube newtonian gravity
//----------------------------------------------------------------------------
{
	double rCube = pow(norm(pos), 3.5);
	return -G/rCube*pos;
}
//============================================================================

//============================================================================
inline vec newtonBeta2(vec pos, vec vel)
//----------------------------------------------------------------------------
// Calculate inverse-square -> inverse-cube newtonian gravity
//----------------------------------------------------------------------------
{
	double rCube = pow(norm(pos), 4);
	return -G/rCube*pos;
}
//============================================================================


//============================================================================
//-------------------------------- MAIN --------------------------------------
//============================================================================
int main(int argc, char *argv[])
//----------------------------------------------------------------------------
// Find the escape velocity of Earth and inverse-square -> inverse-cube gravity
//----------------------------------------------------------------------------
{

	double T = 3;             // simulation time [years]
	int N = 10000;            // N - number of integration points
	int sampleN = 10;         // n - sample point interval

	// Initialize celestial bodies in the Solar System
	Planet Earth(vec({1,0,0}), vec({0,sqrt(2)*2*M_PI,0}),1./333000);
	vector<Planet> solarsystem = vector<Planet>{Earth};

	// Solve motion - escape velocity
	for(int i = 0; i<6; i++)
	{
		solarsystem[0].vel(1) = sqrt(1+0.25*i)*2*M_PI;
		Solver solver(solarsystem, G, true);
		solver.solve(2, newton, T, N, sampleN, "./Raw_Data/data" + to_string(i) + ".txt");
	}
	system("python3 plot.py escapeVel");

	// inverse-square -> inverse-cube gravity
	Solver solverVerlet(solarsystem, G, true);
	solverVerlet.solve(2, newtonBeta1, T, N, sampleN, "./Raw_Data/data.txt");
	system("python3 plot.py beta beta1");

	solverVerlet.solve(2, newtonBeta2, T, N, sampleN, "./Raw_Data/data.txt");
	system("python3 plot.py beta beta2");


	return 0;
}
//============================================================================
