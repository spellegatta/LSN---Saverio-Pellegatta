/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __System__
#define __System__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h> //exit
#include "particle.h"
#include "random.h"

using namespace std;
using namespace arma;

class System {

private:
  const int _ndim = 1;  // Dimensionality of the system
  bool _restart;        // Flag indicating if the simulation is restarted
  int _nblocks;         // Number of blocks for block averaging
  int _nsteps;          // Number of simulation steps in each block
  unsigned int _eqsteps;// Number of steps to thermalize the system
  double _mu, _sigma;   // Mean and dev std of the wave function
  int _nattempts;       // Number of attempted moves
  int _naccepted;       // Number of accepted moves
  double _delta;        // Displacement step
  Random _rnd;          // Random number generator

  // Properties
  int _nprop; // Number of properties being measured
  bool _measure_H_av,_measure_penergy, _measure_kenergy       ;// Flags for measuring average energy, potential energy and kinetic energy
  int _index_H_av, _index_penergy, _index_kenergy;       // Indices for accessing average energy, potential energy and kinetic energy in vec _measurement
  vec _block_av;         // Block averages of properties
  vec _global_av;        // Global averages of properties
  vec _global_av2;       // Squared global averages of properties
  vec _average;          // Average values of properties
  vec _measurement;      // Measured values of properties
  Particle _particle;    // Class representing one particle

public: // Function declarations
  int get_nbl();              // Get the number of blocks
  int get_nsteps();           // Get the number of steps in each block
  void initialize();          // Initialize system properties
  void initialize_properties();// Initialize properties for measurement
  void finalize();            // Finalize system and clean up
  void write_configuration(); // Write final system configuration to XYZ file
    void read_configuration();  // Read system configuration from file
  void write_XYZ(int nconf);  // Write system configuration in XYZ format on the fly
  void initialize_velocities();// Initialize particle velocities
  void step();                // Perform a simulation step
  void block_reset(int blk);  // Reset block averages
  void measure();             // Measure properties of the system
  void averages(int blk);     // Compute averages of properties
  double error(double acc, double acc2, int blk); // Compute error
  void move();        // Move a particle
  void equilibrate(); //Termalizzo il sistema
  bool metro(double, double);       // Perform Metropolis acceptance-rejection step
  void autotuning(); // Aggiusta il delta
};

#endif // __System__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
