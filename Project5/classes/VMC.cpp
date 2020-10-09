//=================================
// Included dependencies
#include "VMC.hpp"
#include <armadillo>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>

using namespace std;
using namespace chrono;
using namespace arma;

//==============
// CONSTRUCTOR
//==============
VMC::VMC(int numDim, int numParticles, myfunc1 acceptAmp, myfunc2 localKinetic,
         myfunc2 localPotential) {
  this->numDim = numDim;
  this->numParticles = numParticles;
  this->acceptAmp = acceptAmp;
  this->localKinetic = localKinetic;
  this->localPotential = localPotential;

  pos = zeros(numDim, numParticles);
  delta = zeros(numDim);
}

//=====================
// Public methods
//=====================

//========================================================
void VMC::mcCycle()
//--------------------------------------------------------
// Brute-force VMC sampling with Metropolis algorithm
//--------------------------------------------------------
{
  for (int j = 0; j < numParticles; j++) {
    for (int k = 0; k < numDim; k++) {
      delta(k) = (myRandu(engine) - 0.5) * step;
    }
    if (myRandu(engine) < acceptAmp(params, omega, pos, delta, j)) {
      pos.col(j) += delta;
      accepted++;
    }
  }
}
//========================================================

//========================================================
Result VMC::solve(int numCycles, int preCycles, double *params, double omega,
                  bool writeToFile)
//--------------------------------------------------------
// Monte Carlo simulation. Update results to struct Result
// and write to file
//--------------------------------------------------------
{
  this->numCycles = numCycles;
  this->preCycles = preCycles;
  this->params = params;
  this->omega = omega;

  ofstream myfile;
  step = 1;
  kineticE = potentialE = k_E = p_E = E = E2 = Var = R12 = 0;
  accepted = 0;

  int interval = int(0.1 * preCycles) + 1;
  for (int i = 0; i < preCycles; i++) {
    mcCycle();

    if ((i + 1) % interval == 0) {
      // tweaks the step if accept rate is more or less than 0.5
      // cout << double(accepted)/(numParticles*interval) << endl;
      step *= 2 * double(accepted) / (numParticles * interval);
      accepted = 0;
    }
  }

  myfile.open("results/data.dat");
  for (int i = 0; i < numCycles; i++) {
    mcCycle();
    // start accumulating data
    kineticE = localKinetic(params, omega, pos);
    potentialE = localPotential(params, omega, pos);
    if (writeToFile) {
      myfile << pos(0, 0) << " " << pos(1, 0) << " " << pos(2, 0) << " "
             << pos(0, 1) << " " << pos(1, 1) << " " << pos(2, 1) << " "
             << kineticE + potentialE << "\n";
    }

    k_E += kineticE;
    p_E += potentialE;
    E += kineticE + potentialE;
    E2 += (kineticE + potentialE) * (kineticE + potentialE);
    R12 += norm(pos.col(0) - pos.col(1));
  }
  myfile.close();

  E /= numCycles;
  E2 /= numCycles;
  Var = E2 - E * E;
  k_E /= numCycles;
  p_E /= numCycles;
  R12 /= numCycles;

  Result myResult = {E, Var, k_E, p_E, R12, accepted};
  return myResult;
}
//========================================================

//========================================================
void VMC::optimize(double *params, int numParams, double range, int step,
                   int maxIter, int numCycles, int preCycles)
//--------------------------------------------------------
// Optimize variational parameters
//--------------------------------------------------------
{
  int paramToChange = numParams - 1; // parameter to increment, start with last
  double energy = 1e10; // minimum energy, starts as an arbitrary large number
  double newEnergy = 0;
  int count = 0; // number of times the parameters have been tweaked
  double *tempParams = new double[2]; // temporary values for the parameters

  while (count < maxIter) {
    tempParams[0] = params[0];
    tempParams[1] = params[1];
    tempParams[paramToChange] -= range;

    for (int i = -step; i <= step; i++) {
      // calculate energy for set of parameters
      newEnergy = this->solve(numCycles, preCycles, tempParams, 1, false).E;
      if (newEnergy < energy) // check if new energy is new minimum
      {
        // set new minimum and parameters
        energy = newEnergy;
        params[paramToChange] = tempParams[paramToChange];
      }
      tempParams[paramToChange] += range / step; // increment parameter
    }

    if (paramToChange > 0) {
      paramToChange -= 1; // change parameter to be incremented
    } else {
      paramToChange = numParams - 1; // swap back to the last parameter
      count++;
      range /= 2; // and decrese the range to increse resolution
      if ((100 * count) % maxIter == 0) {
        cout << (100 * count) / maxIter << '%' << endl;
      }
    }
  }
  delete[] tempParams;
}
//========================================================
