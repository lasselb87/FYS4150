//=================================
// Include guard
#ifndef ISING_HPP
#define ISING_HPP

//=================================
// Included dependencies
#include <cstdlib>
#include <armadillo>
#include <fstream>
#include <cmath>
#include <map>
#include <random>

using namespace std;
using namespace arma;

//=============================================================================
//-------------------------------- CLASS --------------------------------------
//=============================================================================
class Ising
//----------------------------------------------------------------------------
// Initialize a Ising model with periodic boundary conditions
//----------------------------------------------------------------------------
{
private:
	//==========================
	// Declare private variables
	//==========================
	Mat<int> ensemble;               // Matrix of the spins
	int x, y;                        // Coordinates of chosen spin to be flipped
	map<double, double> acceptAmp;   // Map that relates change of energy to
	                                 // the appropirate acceptance amplitude

	//=====================
	// Private methods
	//=====================

	// Initialize parameters
	void init(int L, double T, double J, mt19937_64 &engine);

public:
	//==========================
	// Declare public variables
	//==========================
	int L;                // Lattice Dimenstionality LxL
	double J;             // Coupling constant
	double T;             // Temperature
	double energy;        // Energy of the system
	double magnetization; // Magnetization of the system
	double deltaE;        // Change in energy when flipping single spin
	double deltaM;        // Change in magnetization when flipping single spin

	uniform_int_distribution<int> rand_spin;
	uniform_int_distribution<int> rand_coord;

	//================
	// CONSTRUCTOR(S)
	//================
	Ising();
	Ising(int L, double T, double J, mt19937_64 &engine);
	Ising(Mat<int> ensemble, int L, double T, double J, mt19937_64 &engine);

	//================
	// DESTRUCTOR
	//================
	~Ising();


	//=====================
	// Public methods
	//=====================

	// Periodic boundary conditions
	int periodic(int x);

	// Calculate Ising model energy
	void calcEnergy();

	// Calculate Ising model magnetization
	void calcMagnetization();

	// Print spin matrix
	void print();

	// Trial state for the Metropolis algorithm
	void tryflip(double &aA, mt19937_64 &engine);

	// Flip spin if trial state passes test
	void flip();
};
//=============================================================================

#endif // __ISING_HPP_INCLUDED__
