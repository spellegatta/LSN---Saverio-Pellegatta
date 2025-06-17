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

void write_spin_data(vector<double> M, vector<double> M_error,  vector<double> Cv, vector<double> Cv_error, vector<double> U, vector<double> U_error, vector<double> , vector<double> chi_error, const char file[]);


int main (int argc, char *argv[]){

  int nconf = 1;
  unsigned int n_samples= 20;
  vector<double> Cv, M, Chi, U;
  vector<double> Cv_error, M_error, Chi_error, U_error;
  System SYS;

  cout << endl << "Initializing system..." << endl;
  SYS.initialize();
  SYS.initialize_properties();
  SYS.block_reset(0);
  cout << endl << "SYS initialized" << endl;

  SYS.set_h(0.);
  cout << endl << "h set to 0" << endl;

  for (unsigned int k=0; k<(n_samples+1); k++){
    cout << "# sample: " << k << setw(12) << "T=" << 0.5+(2.0-0.5)*k/double(n_samples) << endl;
    SYS.set_T(0.5+(2.0-0.5)*k/double(n_samples));
    SYS.equilibrate();

    for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
      for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
        SYS.step();
        SYS.measure();

        if(j%50 == 0){
          //SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
          nconf++;
        }
      }

      SYS.averages(i+1);

      if (i==(SYS.get_nbl()-1)){
        Cv.push_back(SYS.get_cv());
        U.push_back(SYS.get_U());
        Chi.push_back(SYS.get_chi());
        Cv_error.push_back(SYS.get_cv_error());
        U_error.push_back(SYS.get_U_error());
        Chi_error.push_back(SYS.get_chi_error());
      }

      SYS.block_reset(i+1);
    }
    SYS.block_reset(0);
    SYS.finalize();
  }

  SYS.initialize();
  SYS.initialize_properties();
  SYS.block_reset(0);
  cout << endl << "SYS re-initialized" << endl;
  SYS.set_h(0.02);
  cout << endl << "h set to 0.02" << endl;
  for (unsigned int k=0; k<(n_samples+1); k++){
    cout << "# sample: " << k << "  T=" << 0.5+(2.0-0.5)*k/double(n_samples) << endl;
    SYS.set_T(0.5+(2.0-0.5)*k/double(n_samples));
    SYS.equilibrate();

    for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
      for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
        SYS.step();
        SYS.measure();
        if(j%50 == 0){
          //SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
          nconf++;
        }
      }
      SYS.averages(i+1);
      if (i==(SYS.get_nbl()-1)){
        M.push_back(SYS.get_M());
        M_error.push_back(SYS.get_M_error());

      }
      SYS.block_reset(i+1);
    }
    SYS.block_reset(0);
    SYS.finalize();
  }
  
  write_spin_data(M, M_error, Cv, Cv_error, U, U_error, Chi, Chi_error, "../OUTPUT/spin_data_metropolis.txt");
  /*
  for (unsigned int i=0; i<M.size(); i++){
    cout << endl << M.at(i) << setw(12) <<M_error.at(i) << setw(12) << Cv.at(i) << setw(12) << Cv_error.at(i) << setw(12) << U.at(i) << setw(12) <<U_error.at(i) << setw(12) << Chi.at(i) << setw(12) << Chi_error.at(i) << endl;
  }
  */

  return 0;
}

void write_spin_data(vector<double> M, vector<double> M_error,  vector<double> Cv, vector<double> Cv_error, vector<double> U, vector<double> U_error, vector<double> chi, vector<double> chi_error, const char file[]){
  ofstream output;
    output.open(file);
    if(output.fail()){
        cout << endl << "Errore di apertura dello stream di output, il file su cui vuoi trascrivere i dati dello spin è nella stessa cartella dell'eseguibile?" << endl;
        return;
    }
    output << "Cv" << setw(12) << "Cv error" << setw(12) << "U" << setw(12) << "U error" << setw(12) << "chi" << setw(12) << "chi error" << setw(12) << "M" << setw(12) << "M error" << endl;

    for(unsigned int i=0; i<M.size(); i++){
        output <<  Cv.at(i) << setw(12) << Cv_error.at(i) << setw(12) <<  U.at(i) << setw(12) << U_error.at(i) << setw(12) <<  chi.at(i) << setw(12) << chi_error.at(i) << setw(12) << M.at(i) << setw(12) << M_error.at(i) << endl;
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
