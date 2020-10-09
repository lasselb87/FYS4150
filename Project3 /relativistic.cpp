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
double GM = 4*M_PI*M_PI;      // Gravitational constant*1 solar mass
double c = 63239.7263;        // Speed of light [AU/yr]

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
inline vec einstein(vec pos, vec vel)
//----------------------------------------------------------------------------
// Calculate relativistic gravity
//----------------------------------------------------------------------------
{
	double r = norm(pos);
	double l = abs(norm(cross(pos, vel)));
	return -GM/pow(r, 3)*(1 + 3*l*l/(r*r*c*c))*pos;
}
//============================================================================


//============================================================================
//-------------------------------- MAIN --------------------------------------
//============================================================================
int main(int argc, char *argv[])
//----------------------------------------------------------------------------
// Calculate the precession of the perihelion of Mercury
//----------------------------------------------------------------------------
{


	double T = 100;            // simulation time [years]
	int N = 1000000000;          // N - number of integration points

	vec pos_mercury = vec({0.3075, 0, 0}); // [AU]
	vec vel_mercury = vec({0, 12.44, 0}); // [AU/yr]
	double M_mercury = 1.65956463e-7; // [solar mass]

	Planet Mercury(pos_mercury, vel_mercury, M_mercury);
	vector<Planet> solarsystem = vector<Planet>{Mercury};

	// Initialize solver
	Solver solver(solarsystem, GM, true);

	// Solve motion
	solver.solvePerihelion(newton, T, N, "./Raw_Data/angle1.txt");
	solver.solvePerihelion(einstein, T, N, "./Raw_Data/angle2.txt");

	system("python3 plot.py precession");

	return 0;
}
//============================================================================
