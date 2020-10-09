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
double G = 4*M_PI*M_PI;        // Gravitational constant


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
//-------------------------------- MAIN --------------------------------------
//============================================================================
int main(int argc, char const *argv[])
//----------------------------------------------------------------------------
// Simulate motion of the Solar System
//
// Command line arguments:
// T - simulation time [years]
// N - number of integration points
// n - sample point interval
//----------------------------------------------------------------------------
{
	// Check if the user is supplying enough commandline arguments
	if (argc < 4)
	{
		cout << "Bad usage! \n";
		cout << "Please supply:\n";
		cout << "Simulation time T\n";
		cout << "Number of time points N\n";
		cout << "Sample point interval sampleN"<< endl;
		return 1;
	}

	double T = atof(argv[1]);
	int N = atoi(argv[2]);
	int sampleN = atoi(argv[3]);

	// Initialize celestial bodies in the Solar System
	Planet Earth(vec({1,0,0}), vec({0,2*M_PI,0}),1./333000);
	vector<Planet> solarsystem = vector<Planet>{Earth};

	// Initialize solver
	Solver solverEuler(solarsystem, G, true);
	Solver solverVerlet(solarsystem, G, true);

	// Solve motion with T=5, N=10 000
	solverEuler.solve(1, newton, T, N, sampleN, "./Raw_Data/data.txt");
	system("python3 plot.py singlePlanet Euler");

	solverVerlet.solve(2, newton, T, N, sampleN, "./Raw_Data/data.txt");
	system("python3 plot.py singlePlanet Verlet");


	// Calculate energies
	ofstream myfile;
	myfile.open("./Raw_Data/energy.txt");
	for(int i=0; i<N/sampleN; i++)
	{
		myfile << setprecision(8)
		       << solverEuler.energyAllPlanets(i) << " "
		       << solverVerlet.energyAllPlanets(i) << "\n";
	}
	myfile.close();
	system("python3 plot.py energy");

	// Calculate energy and angular momentum fluctation with Euler's method
	myfile.open("./Raw_Data/fluctuation.txt");
	for(int n = 100; n<=1e8; n*=10)
	{
		Solver solver(solarsystem, G, true);
		solver.solve(1, newton, T, n, n/10, "./Raw_Data/data.txt");
		myfile << n << " " << solver.totalEnergyFluctuation() << " " <<
		        solver.angularFluctuation() << "\n";
		cout << n << endl;
	}
	myfile.close();
	system("python3 plot.py fluctuation Euler");

	// Calculate energy and angular momentum fluctation with Verlet's method
	myfile.open("./Raw_Data/fluctuation.txt");
	for(int n = 100; n<=1e8; n*=10)
	{
		Solver solver(solarsystem, G, true);
		solver.solve(2, newton, T, n, n/10, "./Raw_Data/data.txt");
		myfile << n << " " << solver.totalEnergyFluctuation() << " " <<
		        solver.angularFluctuation() << "\n";
		cout << n << endl;
	}
	myfile.close();
	system("python3 plot.py fluctuation Verlet");

	return 0;
}
//============================================================================
