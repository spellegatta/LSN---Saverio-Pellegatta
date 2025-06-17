/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#ifndef __Population__ // Direttiva per prevenire inclusioni multiple dello stesso header file
#define __Population__

#include <armadillo> // Inclusione della libreria Armadillo per operazioni con vettori e matrici
#include <fstream>   // Per operazioni sui file (lettura/scrittura)
#include <iomanip>   // Per manipolare l'output (es. setw)
#include <iostream>  // Per input/output su console (cout, cerr)
#include <string>    // Per l'utilizzo delle stringhe
#include <stdlib.h>  // Per la funzione exit()
#include "journey.h" // Inclusione del file di intestazione della classe Journey
#include "random.h"  // Inclusione del file di intestazione per il generatore di numeri casuali

using namespace std;   // Utilizzo del namespace standard C++
using namespace arma;  // Utilizzo del namespace Armadillo

// Dichiarazione della classe Population
class Population {

private:
  unsigned int _ncities;       // Numero di città in ogni percorso di Journey
  unsigned int _npop;          // Numero di percorsi (individui) nella popolazione
  field<Journey> _journey_pop; // Vettore (Armadillo field) che contiene tutti i percorsi Journey della popolazione
  mat _city_coords; // Matrice Armadillo con le coordinate delle città
  Random _rnd;      // Oggetto generatore di numeri casuali
  double _side_radius; // Lato del quadrato (per sim_type=1) o raggio del cerchio (per sim_type=0)
  double _p;           // Esponente per l'algoritmo dell'operatore di selezione
  unsigned int _ngen;  // Numero di generazioni per l'evoluzione
  unsigned int _nContinents; // Numero di "continenti" (processi MPI)
  unsigned int _rank; // Rank del processo MPI corrente
  unsigned int _sim_type; // Tipo di simulazione (0=cerchio, 1=quadrato, 2=da file)
  unsigned int _npeople, _nmigr; // _npeople: numero di individui che migrano; _nmigr: frequenza di migrazione (ogni N generazioni)
  unsigned int _nblocks, _nsteps; // Parametri per l'analisi a blocchi
  double _best_loss; // Il miglior loss (lunghezza del percorso) trovato nella popolazione
  unsigned int _sim_idx; // Indice della simulazione (per il calcolo delle medie a blocchi)

  // Funzioni accessibili solo all'interno della classe (private)
  // Ordina la popolazione di percorsi in base al "loss" (in ordine decrescente o crescente, a seconda dell'implementazione)
  field<Journey> sort_by_loss(field<Journey> field_journeys);

public:
  // Dichiarazioni delle funzioni pubbliche

  // Inizializza l'oggetto System leggendo i parametri dai file di input
  void initialize(unsigned int nContinents, unsigned int rank);
  // Finalizza la simulazione (es. salva il seed del generatore casuale)
  void finalize();
  // Operatore di selezione: restituisce gli indici dei percorsi da selezionare per la riproduzione
  uvec selection_operator(unsigned int npeople);
  // Esegue una ricerca casuale applicando solo operatori di mutazione con una data probabilità
  void rand_search_mutation(double prob);
  // Operatore di crossover: combina due percorsi "genitori" per crearne di nuovi
  void crossover(uvec idx_parents);
  // Evolve il sistema applicando crossover e mutazioni
  void evolution(double p_c, double p_m1, double p_m2, double p_m3,
                 double p_m4);
  // Stampa l'ordine delle città del miglior percorso nel file `cities_order.dat`
  void print_order();
  // Gestisce l'operazione di migrazione tra i continenti (processi MPI)
  void migration();
  // "Olimpiadi": gestisce la comunicazione e l'elezione del miglior percorso tra tutti i continenti
  void Olympics();
  // Crea la popolazione iniziale di percorsi
  void create_pop();
  // Esegue la simulazione utilizzando il metodo delle medie a blocchi
  void blocks_mean();
  // Calcola l'errore statistico per l'analisi a blocchi
  double error(double acc, double acc2, int blk);
  // Restituisce il miglior loss trovato
  double getBestLoss();
};

#endif // __Population__ // Fine della direttiva per prevenire inclusioni multiple

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
