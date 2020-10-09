//=================================
// Included dependencies
#include <cstdlib>
#include <armadillo>
#include <fstream>
#include <cmath>
#include <map>
#include <random>
#include "ising.hpp"
#include "metropolis.hpp"

using namespace std;
using namespace arma;

//==============
// CONSTRUCTOR
//==============
Metropolis::Metropolis(Ising spins)
{
	this->spins = spins;
	rand_float = uniform_real_distribution<float>(0,1);
}

//================
// DESTRUCTOR
//================
Metropolis::~Metropolis(){
	delete[] energyAndMag;
}


//=====================
// Public methods
//=====================

//========================================================
void Metropolis::solve(int cycles, mt19937_64 &engine)
//--------------------------------------------------------
// Metropolis algorithm
//--------------------------------------------------------
{
	energyAndMag = new int[2*cycles];

	energyAndMag[0] = E = spins.energy;
	energyAndMag[cycles] = M = abs(spins.magnetization);

	accepted = 0;
	for(int i=1; i<cycles; i++)
	{
		//Sweeps over LxL spin matrix
		for(int j=0; j<spins.L*spins.L; j++)
		{
			spins.tryflip(acceptAmp, engine);
			if(rand_float(engine) < acceptAmp)
			{
				accepted++;
				spins.flip();
				E += spins.deltaE;
				M += spins.deltaM;
			}
		}
		energyAndMag[i] = E;
		energyAndMag[i+cycles] = abs(M);
	}
}
//========================================================
