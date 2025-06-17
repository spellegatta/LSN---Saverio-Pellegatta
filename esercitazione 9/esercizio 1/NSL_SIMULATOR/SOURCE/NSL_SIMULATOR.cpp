/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>     // Per le operazioni di input/output
#include "population.h" // Include la definizione della classe Population

using namespace std;    // Utilizza lo spazio dei nomi standard

int main (int argc, char *argv[]){

  int nconf = 1; // Variabile non utilizzata nel codice fornito, potrebbe essere un residuo
  Population POP; // Crea un'istanza della classe Population
  POP.initialize(); // Inizializza la popolazione (legge input, crea città, ecc.)
  // Avvia l'evoluzione della popolazione con probabilità di crossover e mutazione
  POP.evolution(0.5, 0.05, 0.05, 0.05, 0.05);
  POP.print_order(); // Stampa l'ordine delle città del miglior percorso trovato

  POP.finalize(); // Finalizza la simulazione (es. salva il seed del generatore casuale)

  return 0; // Il programma termina con successo
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
