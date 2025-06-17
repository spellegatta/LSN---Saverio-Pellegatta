/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Journey__ // Evita inclusioni multiple del file header
#define __Journey__

#include <armadillo> // Include la libreria Armadillo per l'algebra lineare
#include "random.h"  // Include la definizione della classe Random

using namespace std; // Utilizza lo spazio dei nomi standard
using namespace arma; // Utilizza lo spazio dei nomi di Armadillo

// Definizione della classe Journey
class Journey {

  private:
  uvec _cities_order;       // Vettore Armadillo contenente l'ordine delle città nel percorso
  double _loss;             // La "perdita" o costo totale del percorso (ad esempio, la lunghezza)
  unsigned int _ncities;    // Il numero totale di città nel percorso
  Random _rnd;              // Un oggetto generatore di numeri stocastici (casuali)

  public: // Dichiarazioni delle funzioni pubbliche (interfaccia della classe)
  // Costruttore: inizializza un percorso con un dato numero di città, generatore casuale e coordinate
  Journey(unsigned int ncities, Random rnd, mat coord);
  Journey();                // Costruttore di default
  // Operatore di assegnazione: permette di copiare un oggetto Journey in un altro
  Journey& operator=(const Journey& other);

  void compute_loss(mat coord);     // Calcola la perdita (lunghezza) del percorso
  double getLoss();                 // Restituisce il valore della perdita
  uvec getCityOrder();              // Restituisce il vettore che rappresenta l'ordine delle città
  void setCityOrder(uvec cities_order); // Imposta un nuovo ordine per le città

  void check_validity();            // Controlla se il percorso è valido (es. nessuna città visitata due volte)
  void mutation_operator_shift();   // Operatore di mutazione: sposta un blocco di città
  void mutation_operator_swap();    // Operatore di mutazione: scambia due città nel percorso
  void mutation_operator_swap_blocks(); // Operatore di mutazione: scambia due blocchi di città
  void mutation_operator_flip_block();  // Operatore di mutazione: inverte l'ordine di un blocco di città

};

#endif // __Journey__ // Fine della direttiva di inclusione
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
