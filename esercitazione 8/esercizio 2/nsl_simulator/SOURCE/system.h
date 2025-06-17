/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#ifndef __System__ // Include guard to prevent multiple inclusions of this header file.
#define __System__

#include <iostream> // For input/output operations (e.g., cout, cin).
#include <iomanip> // For output formatting (e.g., setw, setprecision).
#include <fstream> // For file stream operations (e.g., ifstream, ofstream).
#include <string> // For string manipulation.
#include <armadillo> // For linear algebra operations (e.g., vec).
#include <stdlib.h> //exit // For standard library functions, including exit.
#include "particle.h" // Include the header file for the Particle class.
#include "random.h" // Include the header file for the Random number generator class.

using namespace std; // Use the standard namespace.
using namespace arma; // Use the Armadillo namespace.

// Definition of the System class, which manages the simulation.
class System {

private: // Private member variables, accessible only within the class.
  int _ndim = 1;    // Dimensionality of the system (initialized to 1).
  bool _restart;        // Flag indicating if the simulation is restarted.
  int _nblocks;         // Number of blocks for block averaging.
  int _nsteps;          // Number of simulation steps in each block.
  unsigned int _eqsteps;// Number of steps to thermalize the system.
  double _mu_new, _sigma_new, _mu_old, _sigma_old;   // Mean and dev std of the wave function before and after metropolis. These are parameters for the trial wave function.
  int _nattempts;       // Number of attempted Monte Carlo moves.
  int _naccepted;       // Number of accepted Monte Carlo moves.
  double _delta;        // Displacement step size for Monte Carlo moves.
  double _H_old, _H_new;//Values for old and new values for H (Hamiltonian or energy expectation).
  double _H_new_err, _H_old_err;//Values for new and oldvalues for H errors. Errors associated with the Hamiltonian estimates.

  double _temp;         //Temperature value, typically used in Simulated Annealing.
  Random _rnd;          // Random number generator object.

  // Properties related to measurements.
  int _nprop; // Number of properties being measured.
  unsigned int _n_bins_psi; // Number of bins for the probability distribution of psi squared.
  double _bin_size_psi, _appo_psi; // Size of each bin and a temporary variable for psi squared.
  bool _measure_H_av,_measure_penergy, _measure_kenergy, _measure_psi     ;// Flags for measuring average energy, potential energy, kinetic energy, and probability distribution of psi.
  int _index_H_av, _index_penergy, _index_kenergy, _index_psi;        // Indices for accessing average energy, potential energy, kinetic energy, and psi squared in the _measurement vector.
  vec _block_av;        // Block averages of properties.
  vec _global_av;       // Global averages of properties (sum of block averages).
  vec _global_av2;      // Squared global averages of properties (sum of squared block averages, for error calculation).
  vec _average;         // Average values of properties for the current block.
  vec _measurement;     // Measured values of properties for the current step.
  Particle _particle;   // Class representing one particle in the system.

public: // Public member functions, accessible from outside the class.
  System(){;} // Default constructor.
  int get_nbl();              // Get the number of blocks.
  int get_nsteps();           // Get the number of steps in each block.
  void initialize();          // Initialize system properties by reading input files.
  void initialize_properties();// Initialize properties for measurement, including setting up output files.
  void set_rnd(Random rnd);         //Initialize _rnd. Set the random number generator object.
  Random get_rnd();                 //Get random generator. Get the current random number generator object.
  void finalize();            // Finalize system and clean up resources.
  void write_configuration(); // Write final system configuration to XYZ file.
  void read_configuration();  // Read system configuration from file.
  void write_XYZ(int nconf);  // Write system configuration in XYZ format on the fly (periodically during simulation).
  void initialize_velocities();// Initialize particle velocities (This function is declared but not defined in the provided .cpp files, suggesting it might be a remnant or intended for future use).
  void step();                // Perform a single simulation step (Monte Carlo move).
  void block_reset(int blk);  // Reset block averages and counters for a new block.
  void measure();             // Measure properties of the system at the current state.
  void averages(int blk);     // Compute and write averages of properties for the current block.
  double error(double acc, double acc2, int blk); // Compute statistical error using block averaging.
  void move();        // Move a particle (propose a Monte Carlo move).
  void equilibrate(bool do_tuning); //Termalizzo il sistema. Equilibrate the system, with an option for autotuning of delta.
  bool metro(double, double);       // Perform Metropolis acceptance-rejection step. Evaluates if a move is accepted based on the Metropolis criterion.
  void autotuning(); // Aggiusta il delta. Adjust the displacement step _delta based on acceptance rate.
  double Boltzmann(); // Calculate Boltzmann factor for Metropolis acceptance (used in Simulated Annealing).
  void SA_proposal(double Rescale_T);     //Calculate new mu and sigma values and make a stime of H. Propose new (mu, sigma) parameters in Simulated Annealing.
  void run_internal_cycle();          //Run the internal metropolis algorithm to evaluate properties. Executes a full simulation cycle (blocks and steps) for a given set of parameters.
  double getT();        //Returns the temperature of the system.
  void setT(double T);        //Set the temperature of the system.
  System(const System &other); // Copy constructor for the System class.
  System& operator=(const System &other); // Assignment operator for the System class.


};

#endif // __System__ // End of include guard.

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
