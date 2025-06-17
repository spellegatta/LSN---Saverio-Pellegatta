/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita" degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <iostream>
#include "system.h" // Include the header file for the System class.

using namespace std;

// Main function where the simulation execution begins.
int main (int argc, char *argv[]){

  int nconf = 1; // Unused variable (possibly for writing configurations, but not implemented in this main).
  double T_limit=1e-6; // Define the lower temperature limit for Simulated Annealing.
  double rescale_T=0.985; // Factor by which the temperature is rescaled in each SA step.
  System SYS; // Create an instance of the System class.
  cout << endl << "Initializing system..." << endl; // Output message to console.
  SYS.initialize(); // Call the initialize method to set up the simulation parameters.
  SYS.initialize_properties(); // Call the initialize_properties method to set up measurement parameters.
  SYS.block_reset(0); // Reset block accumulators before starting.
  cout << endl << "Equilibrating system..." << endl; // Output message to console.
  SYS.equilibrate(true); // Perform the equilibration phase with autotuning enabled.
  cout << endl << "Starting simulation..." << endl; // Output message to console.
  SYS.run_internal_cycle(); // Run the initial simulation cycle to get an estimate for H.
  double T=SYS.getT(); // Get the initial temperature of the system.
  // Start the Simulated Annealing loop.
  while(T>T_limit){ // Continue as long as the temperature is above the limit.
    SYS.SA_proposal(rescale_T); // Propose new parameters and update temperature using Simulated Annealing.
    T=SYS.getT(); // Update the current temperature.
  }

  SYS.finalize(); // Perform finalization steps (e.g., writing final configurations, saving RNG seed).
  cout << endl << "Simulation ended :)" << endl; // Output message indicating simulation completion.

  return 0; // Return 0 to indicate successful execution.
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita" degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
