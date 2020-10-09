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
int main(int argc, char *argv[])
//----------------------------------------------------------------------------
// Simulate the motion of all planets of the Solar System
//----------------------------------------------------------------------------
{


	double T = 250;            // simulation time [years]
	int N = 1000000;          // N - number of integration points
	int sampleN = 10;         // n - sample point interval

	// Initialize celestial bodies in the Solar System
	double M_sun = 1.0;              // [solar mass]
	double M_mercury = 0.166E-6;     // [solar mass]
	double M_venus = 2.081E-6;       // [solar mass]
	double M_earth = 3.039E-6;       // [solar mass]
	double M_mars = 0.323E-6;        // [solar mass]
	double M_jupiter = 954.7E-6;     // [solar mass]
	double M_saturn = 285.8E-6;      // [solar mass]
	double M_uranus = 43.66E-6;      // [solar mass]
	double M_neptune = 51.51E-6;     // [solar mass]
	double M_pluto = 0.006E-6;       // [solar mass]
	vec pos_mercury = vec({-1.323E-01, -4.393E-01, -2.444E-02}); // [AU]
	vec vel_mercury = vec({2.132E-02, 6.577E-03, -2.494E-03});   // [AU/day]
	vel_mercury *= 365.25;                                       // [AU/year]
	vec pos_venus = vec({7.017E-01, 1.868E-01, -3.811E-02});     // [AU]
	vec vel_venus = vec({-5.092E-03, 1.950E-02, 5.612E-04});     // [AU/day]
	vel_venus *= 365.25;                                         // [AU/year]
	vec pos_earth = vec({9.288E-01, 3.700E-01, -9.347E-05});     // [AU]
	vec vel_earth = vec({-6.552E-03, 1.596E-02, 2.174E-08});     // [AU/day]
	vel_earth *= 365.25;                                         // [AU/year]
	vec pos_mars = vec({1.379E+00, -1.325E-01, -3.685E-02});     // [AU]
	vec vel_mars = vec({1.937E-03, 1.512E-02, 2.692E-04});       // [AU/day]
	vel_mars *= 365.25;                                          // [AU/year]
	vec pos_jupiter = vec({-2.654E+00, -4.662E+00, 7.870E-02});  // [AU]
	vec vel_jupiter = vec({6.467E-03, -3.372E-03, -1.306E-04});  // [AU/day]
	vel_jupiter *= 365.25;                                       // [AU/year]
	vec pos_saturn = vec({1.554E+00, -9.934E+00, 1.108E-01});    // [AU]
	vec vel_saturn = vec({5.203E-03, 8.449E-04, -2.219E-04});    // [AU/day]
	vel_saturn *= 365.25;                                        // [AU/year]
	vec pos_uranus = vec({1.717E+01, 1.000E+01, -1.853E-01});    // [AU]
	vec vel_uranus = vec({-2.008E-03, 3.215E-03, 3.792E-05});    // [AU/day]
	vel_uranus *= 365.25;                                        // [AU/year]
	vec pos_neptune = vec({2.892E+01, -7.719E+00, -5.075E-01});  // [AU]
	vec vel_neptune = vec({7.890E-04, 3.051E-03, -8.105E-05});   // [AU/day]
	vel_neptune *= 365.25;                                       // [AU/year]
	vec pos_pluto = vec({1.164E+01, -3.157E+01, 9.700E-03});     // [AU]
	vec vel_pluto = vec({3.008E-03, 4.142E-04, -9.250E-04});     // [AU/day]
	vel_pluto *= 365.25;                                         // [AU/year]

	Planet Sun(vec({0,0,0}), vec({0,0,0}), M_sun);
	Planet Mercury(pos_mercury, vel_mercury, M_mercury);
	Planet Venus(pos_venus, vel_venus, M_venus);
	Planet Earth(pos_earth, vel_earth, M_earth);
	Planet Mars(pos_mars, vel_mars, M_mars);
	Planet Jupiter(pos_jupiter, vel_jupiter, M_jupiter);
	Planet Saturn(pos_saturn, vel_saturn, M_saturn);
	Planet Uranus(pos_uranus, vel_uranus, M_uranus);
	Planet Neptune(pos_neptune, vel_neptune, M_neptune);
	Planet Pluto(pos_pluto, vel_pluto, M_pluto);
	vector<Planet> solarsystem = vector<Planet>{Sun, Mercury, Venus, Earth, Mars,
		                                    Jupiter, Saturn, Uranus, Neptune, Pluto};

	// Set motion properties of the Sun
	for (unsigned int i=1; i<solarsystem.size(); i++) {
		solarsystem[0].pos -= 1/solarsystem[0].M*solarsystem[i].M*solarsystem[i].pos;
		solarsystem[0].vel -= 1/solarsystem[0].M*solarsystem[i].M*solarsystem[i].vel;
	}

	// Initialize solver(s)
	Solver solverVerlet(solarsystem, GM, false);

	// Solve
	solverVerlet.solve(2, newton, T, N, sampleN, "./Raw_Data/data.txt");
	system("python3 plot.py all");


	return 0;
}
//============================================================================
