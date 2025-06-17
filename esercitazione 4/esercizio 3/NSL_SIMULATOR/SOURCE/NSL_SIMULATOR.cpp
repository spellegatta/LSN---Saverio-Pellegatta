/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <iostream>
#include <cmath>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  int nconf = 1;            // counter for configuration snapshots
  System SYS;               // instantiate simulation system

  SYS.initialize();         // read input, set RNG, initial configuration
  SYS.initialize_properties(); // set which properties to measure
  SYS.block_reset(0);       // reset accumulators before data collection

  vec U(SYS.get_nsteps());  // container for internal energy over time
  vec T(SYS.get_nsteps());  // container for temperature over time

  // loop over blocks
  for(int i = 0; i < SYS.get_nbl(); i++){
    SYS.initialize();       // reinitialize system at start of each block
    SYS.initialize_properties();
    SYS.block_reset(0);

    // first half of the block in forward time
    for(int j = 0; j < int(SYS.get_nsteps()/2.); j++){
      SYS.step();           // perform MD/MC step
      SYS.measure();        // measure properties

      if(j % 50 == 0){
        // SYS.write_XYZ(nconf); // save snapshot (disabled)
        nconf++;
      }
      if(i == 0){           // record U and T only in first block
        U.at(j) = SYS.get_U();
        T.at(j) = SYS.get_T();
      }
    }

    SYS.invert_time();      // reverse the velocities/time evolution

    // second half of the block in reversed time
    for(int j = 0; j < ceil(SYS.get_nsteps()/2.); j++){
      SYS.step();           // perform MD/MC step
      SYS.measure();        // measure properties

      if(j % 50 == 0){
        // SYS.write_XYZ(nconf); // save snapshot (disabled)
        nconf++;
      }
      if(i == 0){           // record U and T only in first block
        U.at(j + SYS.get_nsteps()/2) = SYS.get_U();
        T.at(j + SYS.get_nsteps()/2) = SYS.get_T();
      }
    }

    SYS.averages(i + 1);    // compute block averages and errors
    SYS.block_reset(i + 1); // reset before next block
  }

  SYS.write_U_and_T(U, T, 0); // output U and T vs time to file
  SYS.finalize();            // finalize simulation and save seed

  return 0;
}


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
