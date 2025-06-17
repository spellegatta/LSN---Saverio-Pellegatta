/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
// main.cpp

#include <iostream>
#include "system.h"
#include <vector>

using namespace std;

// Write P and U data with their errors to a text file
void write_data(vector<double> P, vector<double> P_error,
                vector<double> U, vector<double> U_error,
                const char file[]);

int main (int argc, char *argv[]){

  int nconf = 1;                         // Counter for configuration snapshots
  unsigned int n_samples = 20;           // Number of temperature samples
  vector<double> P, U;                   // Pressure and internal energy results
  vector<double> P_error, U_error;       // Corresponding errors
  System SYS;                            // Simulation system

  cout << endl << "Initializing system..." << endl;
  SYS.initialize();                      // Read input, initialize RNG and config
  SYS.initialize_properties();           // Set which properties to measure
  SYS.block_reset(0);                    // Reset accumulators before first run
  cout << endl << "SYS initialized" << endl;

  // Loop over a range of temperatures
  for (unsigned int k = 0; k <= n_samples; k++){
    double T = 0.5 + (2.0 - 0.5) * k / double(n_samples);
    cout << "# sample: " << k << setw(12) << "T=" << T << endl;
    SYS.set_T(T);                        // Update system temperature
    SYS.equilibrate();                   // Equilibrate at this T

    // Block-averaging simulation at this T
    for(int i = 0; i < SYS.get_nbl(); i++){
      for(int j = 0; j < SYS.get_nsteps(); j++){
        SYS.step();                      // MD/MC step
        SYS.measure();                   // Measure properties

        if(j % 50 == 0){
          // SYS.write_XYZ(nconf);       // Save snapshot (disabled to avoid I/O overload)
          nconf++;
        }
      }

      SYS.averages(i + 1);               // Compute block averages and errors

      // Store final block's P and U with errors
      if (i == SYS.get_nbl() - 1){
        P.push_back(SYS.get_P());
        U.push_back(SYS.get_U());
        U_error.push_back(SYS.get_U_error());
        P_error.push_back(SYS.get_P_error());
      }

      SYS.block_reset(i + 1);            // Reset for next block
    }

    SYS.block_reset(0);                  // Reset accumulators before next T
    SYS.finalize();                      // Finalize and write out final config/seed
  }

  // Write P and U vs T data to output file
  write_data(P, P_error, U, U_error, "../OUTPUT/U_P_data.txt");

  return 0;
}

void write_data(vector<double> P, vector<double> P_error,
                vector<double> U, vector<double> U_error,
                const char file[]){
  ofstream output;
  output.open(file);
  if(output.fail()){
    cout << endl
         << "Errore di apertura dello stream di output, il file di output non è accessibile."
         << endl;
    return;
  }
  // Header line
  output << "P" << setw(12)
         << "P_error" << setw(12)
         << "U" << setw(12)
         << "U_error" << endl;

  // Data lines
  for(unsigned int i = 0; i < U.size(); i++){
    output << P.at(i) << setw(12)
           << P_error.at(i) << setw(12)
           << U.at(i) << setw(12)
           << U_error.at(i) << endl;
  }
  output.close();
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
