/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#ifndef __Journey__ // Direttiva per evitare inclusioni multiple del file header
#define __Journey__

#include <armadillo> // Include la libreria Armadillo per operazioni con vettori e matrici
#include "random.h"  // Include il file header per il generatore di numeri casuali

using namespace std; // Utilizza il namespace standard
using namespace arma; // Utilizza il namespace Armadillo

// Dichiarazione della classe Journey
class Journey {

private:
  uvec _cities_order;       // Vettore Armadillo che contiene l'ordine delle città nel percorso
  double _loss;             // Valore di "perdita" (costo) del percorso, solitamente la lunghezza totale
  unsigned int _ncities;    // Numero di città nel percorso
  Random _rnd;              // Oggetto generatore di numeri casuali

public:
  // Costruttori:
  // Inizializza un percorso con un dato numero di città, un generatore casuale e le coordinate delle città
  Journey(unsigned int ncities, Random rnd, mat coord);
  // Costruttore di default, inizializza un percorso vuoto
  Journey();
  // Operatore di assegnazione per copiare un oggetto Journey
  Journey &operator=(const Journey &other);

  // Metodi per la gestione del percorso:
  void compute_loss(mat coord); // Calcola la perdita (lunghezza) del percorso basandosi sulle coordinate
  double getLoss();             // Restituisce il valore di perdita del percorso
  uvec getCityOrder();          // Restituisce il vettore con l'ordine delle città
  void setCityOrder(uvec cities_order); // Imposta un nuovo ordine per le città

  // Metodi di validazione e mutazione:
  void check_validity(); // Controlla la validità del percorso (es. ogni città visitata una sola volta)
  void mutation_operator_shift(); // Operatore di mutazione: sposta un blocco di città
  void mutation_operator_swap();  // Operatore di mutazione: scambia due città
  void mutation_operator_swap_blocks(); // Operatore di mutazione: scambia due blocchi di città
  void mutation_operator_flip_block(); // Operatore di mutazione: inverte l'ordine di un blocco di città
};

#endif // Chiusura della direttiva di inclusione

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
