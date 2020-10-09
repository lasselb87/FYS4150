#define CATCH_CONFIG_MAIN

//=================================
// Included dependencies
#include "VMC.hpp"
#include "catch.hpp"
#include "wavefunctions.cpp"
#include <armadillo>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <random>

using namespace std;

//============================================================================
//-------------------------------- TEST --------------------------------------
//============================================================================

double analyticalE(double alpha) {
  return 3. / (2 * alpha) * (1 - alpha * alpha) + 3 * alpha;
}

//=============================================================================
TEST_CASE("Non-interacting benchmark")
//----------------------------------------------------------------------------
// Check that the analytical value for non-interacting energy matches numerical
// for various alpha
//----------------------------------------------------------------------------
{
  double eps = 1e-2;
  VMC solver(3, 2, &acceptAmp1, &localKinetic1, &localPotential1);
  double *params = new double[2];

  for (int i = 0; i <= 10; i++) {
    params[0] = 0.5 + 1. / 10 * i;
    Result myResult = solver.solve(1e6, 1e3, params, 1, false);
    REQUIRE(abs(myResult.E - analyticalE(params[0])) < eps);
  }
}
//=============================================================================

//=============================================================================
TEST_CASE("Acceptence ratio")
//----------------------------------------------------------------------------
// Check that the acceptence ratio is 0.5 within 5%, both with and without
// Jastrow factor
//----------------------------------------------------------------------------
{
  double eps = 2e-2;
  int numCycles = 1e6;
  int preCycles = 1e3;
  int numParticles = 2;
  VMC solver1(3, numParticles, &acceptAmp1, &localKinetic1, &localPotential1);
  VMC solver2(3, numParticles, &acceptAmp2, &localKinetic2, &localPotential2);
  double *params = new double[2];
  params[0] = 1;
  params[1] = 1;

  Result myResult1 = solver1.solve(numCycles, preCycles, params, 1, false);
  Result myResult2 = solver2.solve(numCycles, preCycles, params, 1, false);
  REQUIRE(abs(myResult1.accepted / double(numParticles * numCycles) - 0.5) <
          eps);
  REQUIRE(abs(myResult2.accepted / double(numParticles * numCycles) - 0.5) <
          eps);
}
//=============================================================================
