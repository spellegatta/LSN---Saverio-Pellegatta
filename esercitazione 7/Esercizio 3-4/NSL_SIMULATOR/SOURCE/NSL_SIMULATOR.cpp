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

// Scrive su file i vettori di pressione e energia interna con i rispettivi errori
void write_data(vector<double> P, vector<double> P_error,
                vector<double> U, vector<double> U_error,
                const char file[]);

int main (int argc, char *argv[]){

  int nconf = 1;                         // contatore per i file di configurazione
  vector<double> P, U;                   // vettori per pressione e energia interna medie
  vector<double> P_error, U_error;       // vettori per errori su P e U
  System SYS;                            // istanza dell’oggetto System

  cout << endl << "Initializing system..." << endl;
  SYS.initialize();                      // legge input, inizializza RNG e configurazione
  SYS.initialize_properties();           // setta quali proprietà verranno misurate
  SYS.block_reset(0);                    // resetta accumulatori prima dell’equilibrio
  cout << endl << "SYS initialized" << endl;
  
  cout << endl << "Equilibration of the system..." << endl;
  SYS.equilibrate();                     // fase di equilibrio (warm-up)
  SYS.block_reset(0);                    // reset accumulatore dopo l’equilibrio

  cout << endl << "Starting the simulation..." << endl;
  
  // ciclo sui blocchi di misurazione
  for(int i = 0; i < SYS.get_nbl(); i++){
    // ciclo sui passi (time steps) all’interno di ciascun blocco
    for(int j = 0; j < SYS.get_nsteps(); j++){
      SYS.step();                        // esegue un passo di MD o MC
      SYS.measure();                     // misura le proprietà attuali

      if(j % 50 == 0){
        // SYS.write_XYZ(nconf);         // salva la configurazione in formato XYZ
        nconf++;                         // incremento contatore file
      }
    }

    SYS.averages(i + 1);                 // calcola medie ed errori per il blocco i+1
    SYS.block_reset(i + 1);              // reset accumulatore pronto per il prossimo blocco
  }
  cout << endl << "Simulation completed :)" << endl;
    
  SYS.block_reset(0);                    // reset finale (opzionale)
  SYS.finalize();                        // scrive la configurazione finale e salva il seed

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
