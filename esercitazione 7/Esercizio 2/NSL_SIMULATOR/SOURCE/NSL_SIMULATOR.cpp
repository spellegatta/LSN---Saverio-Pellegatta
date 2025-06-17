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
#include "system.h"
#include <vector>

using namespace std;

// Write pressure and internal energy data with errors to file
void write_data(vector<double> P, vector<double> P_error,
                vector<double> U, vector<double> U_error,
                const char file[]);

int main (int argc, char *argv[]){

  int nconf = 1;          // counter for configuration files
  System SYS;             // instantiate the simulation system

  cout << endl << "Initializing system..." << endl;
  SYS.initialize();       // read input files and set up system
  SYS.initialize_properties(); // enable measurement flags and allocate arrays
  SYS.block_reset(0);     // reset block accumulators before equilibration
  cout << endl << "SYS initialized" << endl;
  
  cout << endl << "Equilibration of the system..." << endl;
  SYS.equilibrate();      // perform equilibration (warm-up) steps

  cout << endl << "Starting the simulation..." << endl;
  // single-block simulation loop (no block averaging here)
  for(int j = 0; j < SYS.get_nsteps(); j++){
    SYS.step();           // perform one MD or MC step
    SYS.measure();        // measure properties at this step

    if(j % 50 == 0){
      // SYS.write_XYZ(nconf); // Write actual configuration in XYZ format
      nconf++;             // increment configuration file counter
    }
  }
  cout << endl << "Simulation completed :)" << endl;
    
  SYS.block_reset(0);     // reset accumulators before finalize
  SYS.finalize();         // write final configuration, save seed, etc.

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
