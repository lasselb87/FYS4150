//=================================
// Include guard
#ifndef METROPOLIS_HPP
#define METROPOLIS_HPP

//=================================
// Included dependencies
#include <cstdlib>
#include <armadillo>
#include <fstream>
#include <cmath>
#include <map>
#include <random>
#include "ising.hpp"

using namespace std;
using namespace arma;

//=============================================================================
//-------------------------------- CLASS --------------------------------------
//=============================================================================
class Metropolis
//----------------------------------------------------------------------------
// Solve the Ising model with the Metropolis algorithm
//----------------------------------------------------------------------------
{
private:
	//==========================
	// Declare private variables
	//==========================
	Ising spins;              // Instance of Ising model
	double acceptAmp;         // Acceptance amplitude

public:
	//==========================
	// Declare public variables
	//==========================
	int* energyAndMag;        // All sampled values of energy and magnetization
	int E;                    // Present energy
	int M;                    // Present magnetization
	int accepted;             // Number of accepted states
	uniform_real_distribution<float> rand_float;  // Random number generator

	//==============
	// CONSTRUCTOR
	//==============
	Metropolis(Ising spins);

	//================
	// DESTRUCTOR
	//================
	~Metropolis();


	//=====================
	// Public methods
	//=====================

	// Metropolis algorithm
	void solve(int cycles, mt19937_64 &engine);
};
//=============================================================================

#endif // __METROPOLIS_HPP_INCLUDED__
