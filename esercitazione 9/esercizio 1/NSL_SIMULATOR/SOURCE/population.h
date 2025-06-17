/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Population__ // Evita inclusioni multiple del file header
#define __Population__

#include <iostream>     // Per input/output standard
#include <iomanip>      // Per manipolare il formato dell'output (es. setw)
#include <fstream>      // Per operazioni su file
#include <string>       // Per l'uso delle stringhe
#include <armadillo>    // Include la libreria Armadillo per l'algebra lineare
#include <stdlib.h>     // Per la funzione exit
#include "journey.h"    // Include la definizione della classe Journey
#include "random.h"     // Include la definizione della classe Random

using namespace std;    // Utilizza lo spazio dei nomi standard
using namespace arma;   // Utilizza lo spazio dei nomi di Armadillo

// Definizione della classe Population (popolazione di percorsi)
class Population {

  private:
    unsigned int _ncities;          // Numero di città in ogni percorso (journey)
    unsigned int _npop;             // Numero di percorsi nella popolazione
    field<Journey> _journey_pop;    // Vettore Armadillo che contiene gli oggetti Journey (i percorsi)
    mat _city_coords;               // Matrice Armadillo con le coordinate X e Y di tutte le città
    Random _rnd;                    // Oggetto generatore di numeri casuali
    double _side_radius;            // Parametro per la disposizione delle città (lato del quadrato o raggio del cerchio)
    double _p;                      // Esponente utilizzato nell'algoritmo dell'operatore di selezione
    unsigned int _ngen;             // Numero di generazioni per l'evoluzione
    unsigned int _sim_type;         // Tipo di simulazione (0 per cerchio, 1 per quadrato per la disposizione delle città)

    // Funzioni accessibili solo all'interno della classe (ausiliarie)
    void create_pop();                                      // Crea la popolazione iniziale di percorsi
    // Ordina la popolazione di percorsi in base alla loro "perdita" (dal peggiore al migliore o viceversa, a seconda dell'implementazione)
    field<Journey> sort_by_loss(field<Journey> field_journeys);

  public: // Dichiarazioni delle funzioni pubbliche (interfaccia della classe)
    void initialize();                                      // Inizializza il sistema leggendo i parametri dai file di input
    void finalize();                                        // Finalizza la simulazione (es. salva il seed)
    // Operatore di selezione: restituisce gli indici dei percorsi che possono generare più discendenti (migliori individui)
    uvec selection_operator(unsigned int npeople);
    // Esegue una ricerca casuale applicando solo operatori di mutazione con una data probabilità
    void rand_search_mutation(double prob);
    void crossover(uvec idx_parents);                       // Operatore di crossover: combina due percorsi parenti per creare nuovi figli
    // Evolve il sistema applicando crossover e diversi tipi di mutazione con le rispettive probabilità
    void evolution(double p_c, double p_m1, double p_m2, double p_m3, double p_m4);
    void print_order();                                     // Stampa l'ordine delle città del miglior percorso trovato

};

#endif // __Population__ // Fine della direttiva di inclusione
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
