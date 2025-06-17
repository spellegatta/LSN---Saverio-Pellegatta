/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita" degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
// main.cpp
#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  int nconf = 1;               // counter for XYZ snapshots
  System SYS;                  // instantiate the system

  cout << endl << "Initializing system..." << endl;
  SYS.initialize();            // read inputs, RNG, initial configuration
  SYS.initialize_properties(); // set up which observables to measure
  SYS.block_reset(0);          // zero out block accumulators

  cout << endl << "Equilibrating system..." << endl;
  SYS.equilibrate();           // perform equilibration phase

  cout << endl << "Starting simulation..." << endl;
  // loop over blocks
  for(int i = 0; i < SYS.get_nbl(); i++){
    // loop over steps within each block
    for(int j = 0; j < SYS.get_nsteps(); j++){
      SYS.step();               // perform one MC step
      SYS.measure();            // measure current observables

      if(j % 50 == 0){
        // SYS.write_XYZ(nconf); // save snapshot (disabled to avoid I/O overload)
        nconf++;
      }
    }
    SYS.averages(i + 1);        // compute and write block averages
    SYS.block_reset(i + 1);     // reset accumulators for next block
  }

  SYS.finalize();               // write final config, save RNG seed
  cout << endl << "Simulation ended :)" << endl;

  return 0;
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
