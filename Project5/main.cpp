//=================================
// Included dependencies
#include "classes/VMC.hpp"
#include "classes/wavefunctions.cpp"
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

//============================================================================
//-------------------------------- MAIN --------------------------------------
//============================================================================
int main(int argc, char const *argv[])
//----------------------------------------------------------------------------
// Calculate variational Monte Carlo estimations with the Metropolis algorithm
//----------------------------------------------------------------------------
{
  double *params1 = new double[2];
  double *params2 = new double[2];
  vector<double> omega;
  //---------------------------------------------------------------------------
  cout << "Benchmarking non-interacting case" << endl;
  VMC solver1(3, 2, &acceptAmp1, &localKinetic1, &localPotential1);
  ofstream myfile;
  myfile.open("results/noninteracting.txt");
  for (int i = 0; i <= 20; i++) {
    params1[0] = 0.5 + 1. / 20 * i;
    Result myResult = solver1.solve(1e6, 1e3, params1, 1, false);
    myfile << params1[0] << " " << myResult.E << " " << myResult.Var << "\n";
  }
  myfile.close();
  //---------------------------------------------------------------------------

  //---------------------------------------------------------------------------
  cout << "Simulation interacting case for different alpha, no Jastrow" << endl;
  VMC solver2(3, 2, &acceptAmp1, &localKinetic1, &localPotential2);
  myfile.open("results/interacting.txt");
  for (int i = 0; i <= 20; i++) {
    params1[0] = 0.5 + 1. / 20 * i;
    Result myResult = solver2.solve(1e6, 1e3, params1, 1, false);
    myfile << params1[0] << " " << myResult.E << " " << myResult.Var << "\n";
  }
  myfile.close();
  //---------------------------------------------------------------------------

  //---------------------------------------------------------------------------
  cout << "Optimizing alpha to get lowest energy for interacting case, no "
          "Jastrow"
       << endl;
  solver2.optimize(params1, 1, 0.6, 10, 10, 1e6, 1e3);
  cout << "Alpha = " << params1[0] << endl;
  //---------------------------------------------------------------------------

  //---------------------------------------------------------------------------
  cout << "Calculating energy for different number of cycles, no Jastrow"
       << endl;
  myfile.open("results/interactingStability.txt");
  for (int i = 1000; i <= 1e6; i *= 10) {
    Result myResult1 = solver2.solve(i, 1e3, params1, 0.05, false);
    Result myResult2 = solver2.solve(i, 1e3, params1, 0.25, false);
    Result myResult3 = solver2.solve(i, 1e3, params1, 1, false);
    myfile << i << " " << myResult1.E << " " << myResult2.E << " "
           << " " << myResult3.E << " " << myResult1.Var << " " << myResult2.Var
           << " "
           << " " << myResult3.Var << "\n";
  }
  myfile.close();
  //---------------------------------------------------------------------------

  //---------------------------------------------------------------------------
  cout << "Calculating mean seperation for different omega" << endl;
  Result myResult1, myResult2, myResult3;
  myfile.open("results/meanSeperation.txt");
  myResult1 = solver2.solve(1e6, 1e3, params1, 0.05, false);
  myResult2 = solver2.solve(1e6, 1e3, params1, 0.25, false);
  myResult3 = solver2.solve(1e6, 1e3, params1, 1, false);
  myfile << myResult1.R12 << " " << myResult2.R12 << " " << myResult3.R12;
  myfile.close();
  //---------------------------------------------------------------------------

  //---------------------------------------------------------------------------
  cout << "Optimizing alpha and beta to get lowest energy, with Jastrow"
       << endl;
  params2[0] = params1[0];
  params2[1] = 1;
  VMC solver3(3, 2, &acceptAmp2, &localKinetic2, &localPotential2);
  solver3.optimize(params2, 2, 0.6, 10, 10, 1e6, 1e3);
  cout << "Alpha = " << params2[0] << endl;
  cout << "Beta = " << params2[1] << endl;
  //---------------------------------------------------------------------------

  //---------------------------------------------------------------------------
  cout << "Calculating mean seperation and energy for different omega, with "
          "Jastrow"
       << endl;
  myfile.open("results/interactionJastrow.txt");
  omega = {0.05, 0.25, 1};
  for (int i = 0; i < omega.size(); i++) {
    Result myResult = solver3.solve(1e7, 1e3, params2, omega[i], false);
    myfile << omega[i] << " " << myResult.E << " " << myResult.Var << " "
           << myResult.R12 << "\n";
  }
  myfile.close();
  //---------------------------------------------------------------------------

  //---------------------------------------------------------------------------
  cout << "Calculating ratio of kinetic and potential energy, both interacting "
          "and non-interacting cast"
       << endl;
  myfile.open("results/virial.txt");
  omega = {0.01, 0.05, 0.1, 0.5, 1};
  params1[0] = 1;
  for (int i = 0; i < omega.size(); i++) {
    myResult1 = solver1.solve(1e6, 1e3, params1, omega[i], false);
    myResult2 = solver3.solve(1e6, 1e3, params2, omega[i], false);
    myfile << omega[i] << " " << myResult1.kinetic / myResult1.potential << " "
           << myResult2.kinetic / myResult2.potential << "\n";
  }
  myfile.close();
  //---------------------------------------------------------------------------

  //---------------------------------------------------------------------------
  cout << "Visualizing the electron cloud" << endl;
  solver1.solve(1e5, 1e3, params1, 1, true);
  //---------------------------------------------------------------------------

  return 0;
}
//============================================================================
