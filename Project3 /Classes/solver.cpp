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
#include "solver.h"
#include "planet.h"
using namespace std;
using namespace arma;

//============================================================================
//------------------------------ METHODS -------------------------------------
//============================================================================


//============================================================================
Solver::Solver(vector<Planet> planets_, double G_, bool staticSun_)
//----------------------------------------------------------------------------
// Constructor
//----------------------------------------------------------------------------
{
	planets = planets_;
	numPlanets = planets.size();
	G = G_;
	staticSun = staticSun_;
}
//============================================================================


//============================================================================
void Solver::totalAcceleration(mat &totalAcc, vec acc(vec, vec))
//----------------------------------------------------------------------------
// Calculate the total acceleration from the force between bodies
//----------------------------------------------------------------------------
{
	totalAcc = zeros(3, numPlanets);
	for(int j=0; j<numPlanets; j++)
	{
		if (staticSun == true) {
			totalAcc.col(j) += acc(pos.col(j), vel.col(j));
		}
		for(int k=j+1; k<numPlanets; k++)
		{
			vec relpos = pos.col(j) - pos.col(k);
			vec relvel = vel.col(j) - vel.col(k);
			mat temp = acc(relpos, relvel);
			totalAcc.col(j) += planets[k].M*temp;
			totalAcc.col(k) -= planets[j].M*temp;
		}
	}
}
//============================================================================


//============================================================================
void Solver::coordinatesToFile(ofstream &myfile)
//----------------------------------------------------------------------------
// Write coordinates to file
//----------------------------------------------------------------------------
{
	for(int j=0; j<numPlanets; j++)
	{
		myfile << pos.col(j)(0) << " "
		       << pos.col(j)(1) << " "
		       << pos.col(j)(2) << " ";
	}
	myfile << "\n";
}
//============================================================================


//============================================================================
void Solver::euler(vec acc(vec, vec))
//----------------------------------------------------------------------------
// Euler's method
//----------------------------------------------------------------------------
{
	prevVel = vel;
	vel = vel + totalAcc*dt;
	pos = pos + prevVel*dt;
	totalAcceleration(totalAcc, acc);
}
//============================================================================


//============================================================================
void Solver::verlet(vec acc(vec, vec))
//----------------------------------------------------------------------------
// Verlet's method
//----------------------------------------------------------------------------
{
	pos = pos + vel*dt + 0.5*totalAcc*dt*dt;
	prevAcc = totalAcc;
	totalAcceleration(totalAcc, acc);
	vel = vel + 0.5*(totalAcc + prevAcc)*dt;
}
//============================================================================


//============================================================================
void Solver::solve(int method, vec acc(vec, vec), double T, int N, int sampleN,
                   string name)
//----------------------------------------------------------------------------
// Solve the system using either Euler's or Verlet's method
//----------------------------------------------------------------------------
{
	solved = true;
	dt = T/N;
	pos = zeros(3, numPlanets);
	vel = zeros(3, numPlanets);

	kineticEnergy      = zeros(N/sampleN, numPlanets);
	potentialEnergy    = zeros(N/sampleN, numPlanets);
	energyAllPlanets   = zeros(N/sampleN);
	angularMomentum    = zeros(N/sampleN);

	for(int i=0; i<numPlanets; i++)
	{
		// Initial conditions
		pos.col(i) = planets[i].pos;
		vel.col(i) = planets[i].vel;
	}

	ofstream myfile;
	myfile.open(name);

	totalAcc = zeros(3, numPlanets);
	prevAcc = zeros(3, numPlanets);
	prevVel = zeros(3, numPlanets);

	totalAcceleration(totalAcc, acc);

	coordinatesToFile(myfile);
	sampleEnergyAndAngular(kineticEnergy, potentialEnergy,
	                       energyAllPlanets, angularMomentum, 0);

	for(int i=0; i<N-1; i++)
	{
		if(method == 1) euler(acc);
		if(method == 2) verlet(acc);
		if ((i+1)%sampleN == 0)
		{
			coordinatesToFile(myfile);
			sampleEnergyAndAngular(kineticEnergy, potentialEnergy,
			                       energyAllPlanets, angularMomentum, (i+1)/sampleN);
		}
	}
	myfile.close();
}
//============================================================================


//============================================================================
void Solver::sampleEnergyAndAngular(mat &kinetic, mat &potential,
                                    vec &energyAllPlanets, vec &angular, int i)
//----------------------------------------------------------------------------
// Calculate energy and angular momentum
//----------------------------------------------------------------------------
{
	vec totalAngular(3,fill::zeros);
	for(int j=0; j<numPlanets; j++)
	{

		kineticEnergy(i,j) = 0.5*planets[j].M*pow(norm(vel.col(j)),2);

		//potential energy from sun
		potentialEnergy(i,j) = -G*planets[j].M/norm(pos.col(j));

		//potential energy inbetween planets
		for(int k=j+1; k<numPlanets; k++)
		{
			double temp = -G*planets[j].M*planets[k].M/
			              norm(pos.col(j) - pos.col(k));
			potentialEnergy(i,j) += temp;
			potentialEnergy(i,k) += temp;
		}
		energyAllPlanets(i) += kineticEnergy(i,j) + potentialEnergy(i,j);
		totalAngular += planets[j].M*cross(pos.col(j), vel.col(j));
	}
	angularMomentum(i) = norm(totalAngular);
}
//============================================================================


//============================================================================
double Solver::kineticFluctuation(int i)
//----------------------------------------------------------------------------
// Find fluctuation of the kinetic energy of planet x
//----------------------------------------------------------------------------
{
	vec kinetic = kineticEnergy.col(i);
	double maxKinetic = kinetic(kinetic.index_max());
	double minKinetic = kinetic(kinetic.index_min());
	return abs((maxKinetic - minKinetic)/maxKinetic);
}
//============================================================================


//============================================================================
double Solver::potentialFluctuation(int i)
//----------------------------------------------------------------------------
// Find fluctuation of the potential energy of planet x
//----------------------------------------------------------------------------
{
	vec potential = potentialEnergy.col(i);
	double maxPotential = potential(potential.index_max());
	double minPotential = potential(potential.index_min());
	return abs((maxPotential - minPotential)/maxPotential);
}
//============================================================================


//============================================================================
double Solver::totalEnergyFluctuation()
//----------------------------------------------------------------------------
// Find fluctuation of the total energy of planet x
//----------------------------------------------------------------------------
{
	double maxEnergy = energyAllPlanets(energyAllPlanets.index_max());
	double minEnergy = energyAllPlanets(energyAllPlanets.index_min());
	return abs((maxEnergy - minEnergy)/maxEnergy);
}
//============================================================================


//============================================================================
double Solver::angularFluctuation()
//----------------------------------------------------------------------------
// Find fluctuation of the angular momentum of planet x
//----------------------------------------------------------------------------
{
	double maxAngular = angularMomentum(angularMomentum.index_max());
	double minAngular = angularMomentum(angularMomentum.index_min());
	return abs((maxAngular - minAngular)/maxAngular);
}
//============================================================================


//============================================================================
void Solver::solvePerihelion(vec acc(vec, vec), double T, int N, string name)
//----------------------------------------------------------------------------
// Calculate perihelion precession
//----------------------------------------------------------------------------
{
	dt = T/N;
	pos = zeros(3, numPlanets);
	vel = zeros(3, numPlanets);

	for(int i=0; i<numPlanets; i++)
	{
		//initial conditions
		pos.col(i) = planets[i].pos;
		vel.col(i) = planets[i].vel;
	}

	ofstream myfile;
	myfile.open(name);

	totalAcc = zeros(3, numPlanets);
	prevAcc = zeros(3, numPlanets);
	vec prevPrevPos = pos.col(0);

	totalAcceleration(totalAcc, acc);
	verlet(acc);
	vec prevPos = pos.col(0);

	verlet(acc);

	for(int i=0; i<N-1; i++)
	{
		if ((norm(prevPos) < norm(prevPrevPos)) and (norm(prevPos) < norm(pos.col(0))))
		{
			myfile << prevPos(1)/prevPos(0) << "\n";
		}
		prevPrevPos = prevPos;
		prevPos = pos.col(0);
		verlet(acc);
	}
	myfile.close();
}
//============================================================================
