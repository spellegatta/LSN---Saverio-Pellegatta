/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <iostream> // Per input/output su console
#include <math.h>   // Per funzioni matematiche, es. sin, cos
#include "journey.h" // Include il file header della classe Journey
#include <cassert>   // Per l'uso della macro assert per debug

using namespace std; // Utilizza il namespace standard

// Costruttore: inizializza un percorso con un numero specificato di città
// e un ordine casuale, partendo sempre dalla città 1.
Journey::Journey(unsigned int ncities, Random rnd, mat coord) {
  _ncities = ncities; // Imposta il numero di città
  _rnd = rnd;         // Assegna il generatore di numeri casuali
  _cities_order.set_size(_ncities); // Alloca la dimensione del vettore per l'ordine delle città

  uvec base = regspace<uvec>(1, _ncities); // Crea un vettore [1, 2, ..., ncities]
  // Estrae tutte le città eccetto la prima (che sarà fissa)
  uvec appo = base.subvec(1, base.size() - 1);
  appo = shuffle(appo); // Mescola casualmente le città rimanenti
  _cities_order.at(0) = base.at(0); // Imposta la prima città (sempre la città 1)
  // Assegna le città mescolate alle posizioni successive
  _cities_order.subvec(1, _cities_order.size() - 1) = appo;

  check_validity();       // Controlla la validità del percorso appena creato
  compute_loss(coord);    // Calcola la lunghezza (perdita) del percorso
}

// Costruttore di default: inizializza un percorso vuoto o non valido
Journey::Journey() {
  _ncities = 0;       // Nessuna città
  Random rnd;         // Crea un generatore di numeri casuali di default
  _rnd = rnd;         // Assegna il generatore
  _cities_order.set_size(_ncities); // Vettore dell'ordine vuoto
  _loss = -1.;        // Perdita non calcolata (segno di stato non valido)
}

// Calcola la lunghezza totale del percorso (la "perdita")
void Journey::compute_loss(mat coord) {
  _loss = 0; // Inizializza la perdita a zero
  // Itera su tutte le città del percorso
  for (unsigned int j = 0; j < coord.n_rows; j++) {
    // Calcola la differenza tra le coordinate della città attuale e della città successiva
    // L'operatore modulo (%) gestisce il ritorno alla prima città per chiudere il percorso
    rowvec diff = coord.row(_cities_order.at(j) - 1) -
                  coord.row((_cities_order.at((j + 1) % _cities_order.size())) - 1);
    _loss += norm(diff, 2); // Aggiunge la distanza euclidea tra le due città alla perdita totale
  }
}

// Restituisce il valore di perdita (lunghezza) del percorso
double Journey::getLoss() { return _loss; }

// Restituisce il vettore che contiene l'ordine delle città nel percorso
uvec Journey::getCityOrder() { return _cities_order; }

// Imposta un nuovo ordine per le città nel percorso e invalida la perdita
void Journey::setCityOrder(uvec cities_order) {
  _cities_order = cities_order; // Assegna il nuovo ordine
  _loss = -1.; // Invalida la perdita, che dovrà essere ricalcolata
  return;
}

// Controlla la validità del percorso:
// 1. Il venditore deve partire dalla città 1.
// 2. Nessuna città deve essere visitata due volte.
// 3. Tutte le città devono essere visitate.
void Journey::check_validity() {
  // Assert: il venditore deve iniziare il suo viaggio dalla città 1.
  assert(("The salesman should start his journey in city 1", _cities_order.at(0) == 1));

  // Controllo per città visitate due volte
  for (unsigned int i = 0; i < _cities_order.size(); i++) {
    for (unsigned int j = i + 1; j < _cities_order.size();
         j++) { // Confronta ogni elemento solo una volta con gli altri
      if (_cities_order.at(i) == _cities_order.at(j)) {
        cerr << endl
             << "A city has been visited twice, exiting..."
             << endl; // Messaggio di errore
        exit(1);     // Termina il programma con codice di errore
      }
    }
  }

  // Controllo per città non visitate
  for (unsigned int i = 0; i < _ncities; i++) {
    // Crea una maschera booleana: 1 se la città (i+1) è presente in _cities_order, 0 altrimenti
    uvec mask = _cities_order == (i + 1);
    // Se 'any(mask)' restituisce false (cioè tutte le città nella maschera sono 0,
    // indicando che la città i+1 non è stata visitata)
    if (!any(mask)) {
      cerr << endl
           << "At least a city has not been visited"
           << endl; // Messaggio di errore
      exit(2);     // Termina il programma con codice di errore
    }
  }

  return;
}

// --- OPERATORI DI MUTAZIONE ---

// Operatore di mutazione: scambia due città casuali nel percorso (escludendo la prima)
void Journey::mutation_operator_swap() {
  unsigned int i, j;
  // Sceglie due indici casuali (escludendo l'indice 0, che corrisponde alla città di partenza fissa)
  i = int(_rnd.Rannyu(1, _ncities));
  j = int(_rnd.Rannyu(1, _ncities));

  // Assicura che gli indici siano positivi (dovrebbero esserlo data la scelta)
  assert(("City order vector indexes should be positive", i > 0 and j > 0));

  // Scambia le città nelle posizioni i e j
  unsigned int appo = _cities_order.at(i);
  _cities_order.at(i) = _cities_order.at(j);
  _cities_order.at(j) = appo;

  check_validity(); // Controlla la validità del percorso dopo la mutazione
  return;
}

// Operatore di mutazione: sposta un blocco di città all'interno del percorso
void Journey::mutation_operator_shift() {
  // 1) Parametri e setup
  const unsigned int N = _cities_order.size(); // Dimensione totale del percorso
  const unsigned int M = N - 1;                // Lunghezza della parte "mobile" (esclusa la prima città)
  uvec pos0 = _cities_order.subvec(0, 0);      // La città fissa di partenza (la prima)
  uvec mobile = _cities_order.subvec(1, N - 1); // La parte del percorso che può essere modificata

  // 2) Estrai parametri casuali per la mutazione
  unsigned int i = unsigned(_rnd.Rannyu(0, M)); // Indice di inizio del blocco da spostare [0..M-1]
  unsigned int m = 1 + unsigned(_rnd.Rannyu(0, M - 1)); // Lunghezza del blocco [1..M]
  // Destinazione del blocco, con "wrap" (se supera la fine, continua dall'inizio)
  unsigned int f = (i + unsigned(_rnd.Rannyu(0, M))) % M;

  // 3) Prendi il blocco di città (anche se "avvolge" la fine del vettore)
  uvec block_vec(m); // Vettore temporaneo per il blocco estratto
  for (unsigned k = 0; k < m; ++k) {
    block_vec(k) = mobile[(i + k) % M]; // Prende le città tenendo conto del "wrap"
  }

  // 4) Prepara il nuovo vettore "mobile" ricostruito
  uvec new_mobile(M);             // Nuovo vettore mobile
  std::vector<bool> filled(M, false); // Vettore booleano per tenere traccia delle posizioni riempite

  // 4a) Inserisci il blocco estratto nella nuova posizione di destinazione, con wrap
  for (unsigned k = 0; k < m; ++k) {
    unsigned pos = (f + k) % M;     // Calcola la posizione con wrap
    new_mobile(pos) = block_vec(k); // Inserisce la città del blocco
    filled[pos] = true;             // Segna la posizione come riempita
  }

  // 4b) Rimetti in ordine tutti gli altri elementi che non facevano parte del blocco estratto
  unsigned rest_idx = (i + m) % M; // Primo elemento "non-bloccato" nella vecchia mobile
  for (unsigned pos = 0; pos < M; ++pos) {
    if (!filled[pos]) { // Se la posizione non è stata riempita dal blocco
      new_mobile(pos) = mobile[rest_idx]; // Inserisce l'elemento rimanente
      rest_idx = (rest_idx + 1) % M;      // Passa al prossimo elemento rimanente
    }
  }

  // 5) Ricostruisci l'ordine completo del percorso
  _cities_order.set_size(N);         // Riassegna la dimensione totale
  _cities_order[0] = pos0(0);        // Reinserisce la città fissa di partenza
  _cities_order.subvec(1, N - 1) = new_mobile; // Inserisce il nuovo vettore mobile

  // 6) Check finale
  check_validity(); // Controlla la validità del percorso dopo la mutazione
}

// Operatore di mutazione: scambia due blocchi di città nel percorso
void Journey::mutation_operator_swap_blocks() {
  // 1) Dimensioni
  unsigned int N = static_cast<unsigned int>(_cities_order.n_elem); // Numero totale di elementi
  unsigned int M = N - 1; // Lunghezza della parte "mobile" (esclusa la prima città)

  // 2) Vista diretta sulla sub-vettore [1..N-1] (la parte mutabile)
  uvec mobile = _cities_order.subvec(1, N - 1);

  // 3) Scegli la lunghezza del blocco in [1..M/2]
  unsigned int block_length = 1 + static_cast<unsigned int>(_rnd.Rannyu(0, M / 2));

  // 4) Scegli l'indice di partenza del primo blocco in [0..M-1]
  unsigned int start1 = static_cast<unsigned int>(_rnd.Rannyu(0, M));

  // 5) Scegli di quanto spostare il secondo blocco rispetto al primo in [1..M-1]
  unsigned int nshifts = 1 + static_cast<unsigned int>(_rnd.Rannyu(0, M - 1));

  // 6) Esegui lo swap ciclico (con modulo M per gestione del wrap-around)
  for (unsigned int k = 0; k < block_length; ++k) {
    unsigned int idx1 = (start1 + k) % M; // Indice della città nel primo blocco
    // Indice della città corrispondente nel secondo blocco (spostato di nshifts)
    unsigned int idx2 = (start1 + nshifts + k) % M;
    swap(mobile(idx1), mobile(idx2)); // Scambia le due città
  }

  // 7) La subvec modifica direttamente _cities_order, quindi basta validare
  check_validity(); // Controlla la validità del percorso dopo la mutazione
}

// Operatore di mutazione: inverte l'ordine di un blocco di città
void Journey::mutation_operator_flip_block() {
  // Sceglie una lunghezza casuale per il blocco da invertire (escludendo la prima città)
  unsigned int block_length = (unsigned int)(_rnd.Rannyu(1, int(_cities_order.size() - 1)));
  uvec mobile = _cities_order.subvec(1, _cities_order.size() - 1); // Parte mutabile del percorso
  uvec pos0 = _cities_order.subvec(0, 0); // La città di partenza fissa

  // Sceglie un indice di inizio casuale per il blocco nella parte mobile
  unsigned int start = (unsigned int)(_rnd.Rannyu(0, mobile.size()));

  // Caso 1: il blocco non "avvolge" la fine del vettore
  if (start + block_length <= mobile.size()) {
    uvec appo = mobile.subvec(start, start + block_length - 1); // Estrae il blocco
    appo = flipud(appo); // Inverte l'ordine del blocco
    mobile.subvec(start, start + block_length - 1) = appo; // Reinserisce il blocco invertito
  } else {
    // Caso 2: il blocco "avvolge" la fine del vettore (es. inizia alla fine e finisce all'inizio)
    uvec appo1 = mobile.subvec(start, mobile.size() - 1); // Prima parte del blocco (fino alla fine)
    // Seconda parte del blocco (dall'inizio del vettore)
    uvec appo2 = mobile.subvec(0, start + block_length - mobile.size() - 1);
    uvec appo = flipud(join_vert(appo1, appo2)); // Unisce le due parti e inverte l'intero blocco

    // Reinserisce le parti invertite nel vettore mobile originale
    mobile.subvec(start, mobile.size() - 1) =
        appo.subvec(0, mobile.size() - 1 - start);
    mobile.subvec(0, start + block_length - mobile.size() - 1) =
        appo.subvec(mobile.size() - start, appo.size() - 1);
  }
  // Ricostruisce l'ordine completo unendo la città fissa e il vettore mobile modificato
  _cities_order = join_vert(pos0, mobile);

  check_validity(); // Controlla la validità del percorso dopo la mutazione
}

// Operatore di assegnazione (copy assignment operator)
Journey &Journey::operator=(const Journey &other) {
  // cout << endl <<"before cities " << other._ncities << " order " <<
  // other._cities_order.t() << " loss " << other._loss << endl; // Debugging

  if (this != &other) { // Evita l'auto-assegnazione
    _ncities = other._ncities;       // Copia il numero di città
    _rnd = other._rnd;               // Copia il generatore di numeri casuali
    _cities_order = other._cities_order; // Copia l'ordine delle città
    _loss = other._loss;             // Copia la perdita
  }
  // cout << endl <<"after cities " << _ncities << " order " <<
  // _cities_order.t() << " loss " << _loss << endl; // Debugging
  return *this; // Restituisce un riferimento all'oggetto corrente
}