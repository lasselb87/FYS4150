//=================================
// Include guard
#ifndef PLANET_H
#define PLANET_H

//=================================
// Included dependencies
#include <armadillo>
using namespace std;
using namespace arma;

//=============================================================================
//-------------------------------- CLASS --------------------------------------
//=============================================================================
class Planet
//----------------------------------------------------------------------------
// Initialize planet with properties 'position', 'velocity', and 'mass'
// Energies and angular momentum is also calculated.
//----------------------------------------------------------------------------
{
public:
	//==========================
	// Declare public variables
	//==========================
	vec pos;           // Position vector
	vec vel;           // Velocity vector
	double M;          // Mass of planet

	//================
	// Constructor(s)
	//================
	Planet();
	Planet(vec position, vec velocity, double mass);
};
//=============================================================================

#endif // __PLANET_H_INCLUDED__
