/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <iostream>    // Per input/output su console
#include "population.h" // Inclusione del file di intestazione della classe Population
#include <mpi.h>       // Per le funzionalità MPI (Message Passing Interface)

using namespace std;   // Utilizzo del namespace standard C++

// Funzione main: punto di ingresso del programma
int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv); // Inizializza l'ambiente MPI
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size); // Ottiene il numero totale di processi MPI (continenti)
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Ottiene il rank (ID) del processo corrente

  // Solo il processo con rank 0 stampa un messaggio di inizializzazione
  if (rank == 0)
    cout << size << " nodes initialized!" << endl;

  unsigned int N_migr = 10; // Variabile non usata in questo main, ma definita
  unsigned int N_people = 10; // Variabile non usata in questo main, ma definita
  int nconf = 1;              // Variabile non usata in questo main, ma definita

  Population POP; // Crea un oggetto Population
  // Inizializza l'oggetto POP con il numero di continenti e il rank del processo corrente
  POP.initialize(size, rank);

  POP.blocks_mean(); // Esegue la simulazione utilizzando il metodo delle medie a blocchi

  POP.finalize(); // Finalizza la simulazione (es. salva il seed del generatore casuale)

  MPI_Finalize(); // Termina l'ambiente MPI
  return 0;       // Il programma termina con successo
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
