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
#include <vector>
#include <stdlib.h> //exit
#include "particle.h"
#include "random.h"

using namespace std;
using namespace arma;

class System {

private:
  const int _ndim = 3;              // Dimensionality of the system
  bool _restart;                    // Flag indicating if the simulation is restarted
  int _sim_type;                    // Type of simulation (e.g., Lennard-Jones, Ising)
  int _npart;                       // Number of particles
  int _nblocks;                     // Number of blocks for block averaging
  int _nsteps;                      // Number of simulation steps in each block
  int _nattempts;                   // Number of attempted moves
  int _naccepted;                   // Number of accepted moves
  double _temp, _beta;              // Temperature and inverse temperature
  double _rho, _volume;             // Density and volume of the system
  double _r_cut;                    // Cutoff radius for pair interactions
  double _delta;                    // Displacement step for particle moves
  double _J, _H;                    // Parameters for the Ising Hamiltonian
  vec _side;                        // Box dimensions
  vec _halfside;                    // Half of box dimensions
  Random _rnd;                      // Random number generator
  field<Particle> _particle;        // Field of particle objects representing the system
  vec _fx, _fy, _fz;                // Forces on particles along x, y, and z directions
  unsigned int _n_part_excluded = 0; // Particles excluded in pofv count

  // Properties
  unsigned int _neq;               // Number of equilibration steps
  int _nprop;                      // Number of properties being measured
  bool _measure_penergy, _measure_kenergy, _measure_tenergy;  // Flags for energy measurements
  bool _measure_temp, _measure_pressure, _measure_gofr;        // Flags for T, P, and g(r)
  bool _measure_magnet, _measure_cv, _measure_chi;             // Flags for M, Cv, χ
  bool _measure_pofv;                                          // Flag for velocity distribution
  int _index_penergy, _index_kenergy, _index_tenergy, _index_E1, _index_E2;  // Indices into measurement vector
  int _index_temp, _index_pressure, _index_gofr;               // Indices for T, P, g(r)
  int _index_magnet, _index_cv, _index_chi;                    // Indices for M, Cv, χ
  int _index_pofv;                                             // Index for velocity distribution
  int _n_bins;           // Number of bins for radial distribution function
  int _n_bins_v;         // Number of bins for velocity distribution
  double _bin_size;      // Bin size for radial distribution function
  double _bin_size_v;    // Bin size for velocity distribution
  double _vtail, _ptail; // Tail corrections for energy and pressure
  double _v_max;         // Maximum velocity for pofv
  vec _block_av;         // Block accumulators of properties
  vec _global_av;        // Global averages (sum of blocks)
  vec _global_av2;       // Global averages of squares
  vec _average;          // Average values per block
  vec _measurement;      // Instantaneous measured values
  vector<double> _cv, _M, _U, _chi, _P; // Results and errors for Cv, M, U, χ, P

public: // Function declarations
  int get_nbl();               // Get the number of blocks
  int get_nsteps();            // Get the number of steps in each block
  void initialize();           // Initialize system properties and configuration
  void initialize_properties(); // Set up which properties to measure
  void finalize();             // Finalize the simulation and write outputs
  void write_configuration();  // Write final system configuration to XYZ
  void write_XYZ(int nconf);   // Write configuration nconf on the fly
  void read_configuration();   // Read configuration from file
  void initialize_velocities(); // Initialize velocities for MD
  void step();                 // Perform a simulation step (MD or MC)
  void block_reset(int blk);   // Reset block averages before a new block
  void measure();              // Measure current properties
  void averages(int blk);      // Compute and write averages for block blk
  void set_h(double h);        // Set external field H for Ising
  void set_T(double T);        // Set temperature T
  void equilibrate();          // Perform equilibration steps
  double error(double acc, double acc2, int blk); // Compute statistical error
  void move(int part);         // Attempt a Monte Carlo move on part
  bool metro(int part);        // Metropolis accept/reject step
  double pbc(double position, int i); // Apply periodic boundary conditions to coordinate
  int pbc(int i);              // Apply PBC to spin index
  void Verlet();               // Perform Verlet integration (MD)
  double Force(int i, int dim); // Compute force on particle i along dim
  double Boltzmann(int i, bool xnew); // Compute Metropolis Boltzmann factor
  double get_cv();             // Return computed heat capacity
  double get_chi();            // Return computed susceptibility
  double get_M();              // Return computed magnetization
  double get_U();              // Return computed total energy
  double get_P();              // Return computed pressure
  double get_cv_error();       // Return error on Cv
  double get_chi_error();      // Return error on χ
  double get_M_error();        // Return error on M
  double get_U_error();        // Return error on U
  double get_P_error();        // Return error on P

};

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
