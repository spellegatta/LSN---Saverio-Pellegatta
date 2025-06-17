/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <cmath>    // Per funzioni matematiche (es. pow, sqrt)
#include <cstdlib>  // Per la funzione exit()
#include <string>   // Per l'utilizzo delle stringhe (es. to_string)
#include "population.h" // Inclusione del file di intestazione della classe Population
#include <cassert>  // Per l'utilizzo della macro assert
#include <mpi.h>    // Per le funzionalità MPI (Message Passing Interface)

using namespace std;  // Utilizzo del namespace standard C++
using namespace arma; // Utilizzo del namespace Armadillo

// Metodo per inizializzare l'oggetto Population e i parametri della simulazione
void Population ::initialize(unsigned int nContinents,
                            unsigned int rank) { // Inizializza l'oggetto System secondo il contenuto dei file di input

  int p1, p2; // Variabili per leggere due numeri primi per l'inizializzazione del RNG
  ifstream Primes("../INPUT/Primes"); // Apre il file Primes
  Primes >> p1 >> p2; // Legge i due numeri primi
  Primes.close();     // Chiude il file

  int seed[4]; // Array per il seed del generatore di numeri casuali
  ifstream Seed("../INPUT/seed.in"); // Apre il file seed.in
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3]; // Legge il seed
  // Imposta il generatore di numeri casuali con il seed e numeri primi unici per ogni rank MPI
  _rnd.SetRandom(seed, p1 + rank, p2 + rank);
  _sim_idx = 0; // Inizializza l'indice della simulazione

  _nContinents = nContinents; // Inizializza il numero di continenti (processi MPI)
  _rank = rank;               // Inizializza il rank del processo MPI corrente

  ifstream input("../INPUT/input.dat"); // Apre il file di input principale
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat"); // Apre il file di output generale

  // Se il rank è 0 (il processo "master"), apre un file per i migliori loss a blocchi
  if (_rank == 0) {
    ofstream coutl("../OUTPUT/mb_best_loss.dat");
    if (coutl.is_open()) {
      coutl << "# BLOCKS: " << setw(12) << "BLOCK_AVERAGE: " << setw(12)
            << "SIMULATION_AVERAGE: " << setw(12)
            << "ERROR: " << endl; // Scrive l'intestazione del file
    } else {
      cerr << endl
           << "Couln't open mb_best_loss.dat"
           << endl; // Errore se il file non può essere aperto
    }
    coutl.close(); // Chiude il file
  }

  string property; // Variabile temporanea per leggere i nomi delle proprietà dal file di input
  double delta;    // Variabile non usata in questo contesto
  // Ciclo per leggere le proprietà dal file di input
  while (!input.eof()) {
    input >> property; // Legge il nome della proprietà
    if (property == "NCITIES") {
      input >> _ncities; // Legge il numero di città
    } else if (property == "SIDE_RADIUS") {
      input >> _side_radius; // Legge il lato/raggio della configurazione
    } else if (property == "SIMULATION_TYPE") {
      input >> _sim_type; // Legge il tipo di simulazione
      if (_sim_type > 2) {
        cerr << "PROBLEM: unknown simulation type"
             << endl; // Errore per tipo di simulazione sconosciuto
        exit(EXIT_FAILURE); // Esce dal programma
      }
      if (_sim_type == 0) { // Città su una circonferenza
        coutf << "CITIES ON A CIRCULAR CONFIGURATION"
              << endl; // Scrive nel file di output
        _city_coords.set_size(_ncities, 2); // Ridimensiona la matrice delle coordinate
        ofstream coutc;
        coutc.open("../OUTPUT/coordinates.dat"); // Apre il file per le coordinate
        coutc << "X:     Y: " << endl;           // Intestazione
        for (unsigned int j = 0; j < _ncities; j++) {
          double theta = _rnd.Rannyu(0, 2 * M_PI); // Genera un angolo casuale
          _city_coords(j, 1) = _side_radius * sin(theta); // Calcola coordinata Y
          _city_coords(j, 0) = _side_radius * cos(theta); // Calcola coordinata X
          coutc << _city_coords(j, 0) << setw(12) << _city_coords(j, 1)
                << endl; // Scrive le coordinate
        }
        coutc.close(); // Chiude il file
      } else if (_sim_type == 1) { // Città in un'area quadrata
        _city_coords.set_size(_ncities, 2);
        coutf << "CITIES IN A SQUARE AREA" << endl;
        ofstream coutc;
        coutc.open("../OUTPUT/coordinates.dat");
        coutc << "X:     Y: " << endl;
        for (unsigned int j = 0; j < _ncities; j++) {
          _city_coords(j, 1) =
              _rnd.Rannyu(0, _side_radius); // Genera coordinata Y casuale
          _city_coords(j, 0) =
              _rnd.Rannyu(0, _side_radius); // Genera coordinata X casuale
          coutc << _city_coords(j, 1) << setw(12) << _city_coords(j, 0)
                << endl;
        }
        coutc.close();
      } else if (_sim_type == 2) { // Città caricate da un file
        _city_coords.set_size(_ncities, 2);
        ifstream citycoords_stream("../INPUT/cap_prov_ita.dat"); // Apre il file delle coordinate
        if (!citycoords_stream.is_open()) {
          cerr << endl
               << "Couldn't open ../INPUT/cap_prov_ita.dat, exiting..."
               << endl; // Errore se il file non può essere aperto
          exit(5);
        } else {
          for (unsigned int i = 0; i < _ncities; i++) {
            citycoords_stream >> _city_coords(i, 0) >>
                _city_coords(i, 1); // Legge le coordinate
          }
          unsigned int appo;
          citycoords_stream >>
              appo; // Tenta di leggere un altro valore per verificare EOF
          if (!citycoords_stream.eof()) {
            cerr << endl
                 << "City coordinates not loaded correctly, exiting..."
                 << endl; // Errore se le coordinate non sono caricate correttamente
            exit(6);
          }
        }
      }
    } else if (property == "NPOP") {
      input >> _npop; // Legge la dimensione della popolazione
    } else if (property == "NBLOCKS") {
      input >> _nblocks; // Legge il numero di blocchi per l'analisi statistica
    } else if (property == "NSTEPS") {
      input >> _nsteps; // Legge il numero di passi per blocco
    } else if (property == "NGEN") {
      input >> _ngen; // Legge il numero di generazioni
    } else if (property == "NMIGR") {
      input >> _nmigr; // Legge la frequenza di migrazione
    } else if (property == "NPEOPLE") {
      input >> _npeople; // Legge il numero di individui che migrano
    } else if (property == "SELECTION_P") {
      input >> _p; // Legge l'esponente per l'operatore di selezione
    } else if (property == "ENDINPUT") {
      coutf << "Reading input completed!" << endl; // Messaggio di completamento
      break;                                      // Esce dal ciclo di lettura
    } else
      cerr << "PROBLEM: unknown input" << endl; // Errore per input sconosciuto
  }
  input.close();                      // Chiude il file di input
  _journey_pop.set_size(_npop);       // Ridimensiona la popolazione di percorsi
  create_pop();                       // Crea la popolazione iniziale
  coutf << "System initialized!" << endl; // Messaggio di inizializzazione completata
  coutf.close();                      // Chiude il file di output generale
  return;
}

// Metodo per gestire la migrazione di individui tra i "continenti" (processi MPI)
void Population ::migration() {
  unsigned int c1, c2; // Indici dei continenti coinvolti nella migrazione

  // Solo il rank 0 (master) sceglie casualmente i due continenti per la migrazione
  if (_rank == 0) {
    c1 = static_cast<unsigned int>(_rnd.Rannyu(0, _nContinents)); // Sceglie il primo continente
    do {
      c2 = static_cast<unsigned int>(_rnd.Rannyu(0, _nContinents)); // Sceglie il secondo continente
    } while (c2 == c1); // Assicura che i due continenti siano diversi
  }

  // Trasmette gli indici dei continenti scelti a tutti i processi
  MPI_Bcast(&c1, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(&c2, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  // Solo i continenti coinvolti nella migrazione procedono
  if (_rank == c1 || _rank == c2) {
    uvec migrants_idxs; // Vettore per gli indici dei migranti

    // Se il rank corrente è il primo continente (c1)
    if (_rank == c1) {
      // Seleziona gli indici degli individui che migreranno
      migrants_idxs = selection_operator(_npeople);
      // Invia gli indici dei migranti al secondo continente (c2)
      MPI_Send(migrants_idxs.memptr(), (int)_npeople, MPI_UNSIGNED_LONG, c2,
               777, MPI_COMM_WORLD);
    }
    // Se il rank corrente è il secondo continente (c2)
    else if (_rank == c2) {
      migrants_idxs.set_size(_npeople); // Ridimensiona il vettore per ricevere gli indici
      // Riceve gli indici dei migranti dal primo continente (c1)
      MPI_Recv(migrants_idxs.memptr(), (int)_npeople, MPI_UNSIGNED_LONG, c1,
               777, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Crea un "field" di percorsi (viaggi) per i migranti
    field<uvec> migrants(migrants_idxs.size());
    // Copia l'ordine delle città dei migranti selezionati
    for (unsigned int i = 0; i < migrants.size(); i++) {
      migrants.at(i) = _journey_pop.at(migrants_idxs.at(i)).getCityOrder();
    }

    int itag1 = 1; // Tag per il messaggio di invio
    int itag2 = 2; // Tag per il messaggio di ricezione
    MPI_Status stat1, stat2;
    MPI_Request req;

    // Ciclo per scambiare i percorsi completi dei migranti
    for (unsigned int i = 0; i < migrants.size(); i++) {
      if (_rank == c1) { // Il primo continente invia e poi riceve
        unsigned int n = migrants.at(i).size();
        // Invia in modo non bloccante il percorso del migrante a c2
        MPI_Isend(migrants.at(i).memptr(), n, MPI_UNSIGNED_LONG, c2, itag1,
                  MPI_COMM_WORLD, &req);
        // Riceve il percorso del migrante da c2
        MPI_Recv(migrants.at(i).memptr(), n, MPI_UNSIGNED_LONG, c2, itag2,
                 MPI_COMM_WORLD, &stat2);
        // Attende che l'invio non bloccante sia completato
        MPI_Wait(&req, &stat1);
      } else if (_rank == c2) { // Il secondo continente riceve e poi invia
        unsigned int n = migrants.at(i).size();
        // Invia il percorso del migrante a c1
        MPI_Send(migrants.at(i).memptr(), n, MPI_UNSIGNED_LONG, c1, itag2,
                 MPI_COMM_WORLD);
        // Riceve il percorso del migrante da c1
        MPI_Recv(migrants.at(i).memptr(), n, MPI_UNSIGNED_LONG, c1, itag1,
                 MPI_COMM_WORLD, &stat1);
      }
    }
    // Aggiorna la popolazione locale con i percorsi dei migranti ricevuti
    for (unsigned int i = 0; i < migrants.size(); i++) {
      _journey_pop.at(migrants_idxs(i)).setCityOrder(migrants.at(i));
      _journey_pop.at(migrants_idxs(i)).compute_loss(_city_coords);
    }
    // Riordina la popolazione dopo la migrazione
    _journey_pop = sort_by_loss(_journey_pop);
    return;
  } else {
    return; // I continenti non coinvolti non fanno nulla
  }
}

// Metodo per creare la popolazione iniziale di percorsi (Journey)
void Population ::create_pop() {
  for (unsigned int i = 0; i < _npop; i++) {
    // Crea un nuovo oggetto Journey con il numero di città, il generatore casuale e le coordinate
    Journey appo_journey(_ncities, _rnd, _city_coords);
    _journey_pop.at(i) = appo_journey; // Aggiunge il percorso alla popolazione
  }
  // Ordina la popolazione in base al loss (dal peggiore al migliore, a seconda di sort_by_loss)
  _journey_pop = sort_by_loss(_journey_pop);

  return;
}

// Metodo per ordinare un "field" di oggetti Journey in base al loro "loss"
field<Journey> Population ::sort_by_loss(field<Journey> field_pop) {
  vec appo_losses(field_pop.size()); // Vettore per memorizzare i valori di loss

  // Estrae i loss da ogni Journey nel field_pop
  for (unsigned int i = 0; i < field_pop.size(); i++) {
    appo_losses.at(i) = field_pop.at(i).getLoss();
  }

  // Ordina gli indici in base ai loss, in ordine "discendente" (che in questo contesto significa dal loss più grande al più piccolo,
  // cioè dal percorso peggiore al migliore)
  uvec idx = sort_index(appo_losses, "descend");
  field<Journey> appo(_npop); // Crea un nuovo field per la popolazione ordinata

  // Ricostruisce il field_pop in base all'ordine degli indici
  for (unsigned int j = 0; j < field_pop.size(); j++) {
    appo.at(j) = field_pop.at(idx.at(j));
  }
  field_pop = appo; // Assegna la popolazione ordinata
  return field_pop; // Restituisce la popolazione ordinata
}

// Metodo per finalizzare la simulazione
void Population ::finalize() {
  _rnd.SaveSeed(); // Salva lo stato del generatore di numeri casuali
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat", ios::app); // Apre il file di output in modalità append
  coutf << "Simulation completed!" << endl;      // Scrive un messaggio di completamento
  coutf.close();                                 // Chiude il file
  return;
}

/////////////////////////////SELECTION OPERATOR///////////////////////////////////////////////////

// Operatore di selezione: restituisce gli indici degli individui selezionati
// per la riproduzione, privilegiando quelli con un loss migliore (più basso)
uvec Population ::selection_operator(unsigned int n_people) {
  uvec idx(n_people); // Vettore per memorizzare gli indici selezionati
  for (unsigned int i = 0; i < n_people; i++) {
    // La selezione è basata su una distribuzione di probabilità che favorisce
    // gli individui con indici più bassi (e quindi loss migliori, se la popolazione è ordinata)
    // Formula: index = NPOP * (random_number)^p
    idx.at(i) = (unsigned int)(_npop * pow(_rnd.Rannyu(), _p));
  }
  return idx; // Restituisce il vettore degli indici selezionati
}

/////////////////////////////////////RANDOM SEARCH//////////////////////////////////////////////////
// Questa funzione esegue una ricerca casuale applicando solo mutazioni
void Population ::rand_search_mutation(double prob) {
  ofstream coutl("../OUTPUT/RS_mutation.dat"); // Apre il file di output per la mutazione RS
  if (coutl.is_open()) {
    coutl << "GEN:       BEST_LOSS:       AVERAGE_LOSS: "
          << endl; // Intestazione del file
  } else {
    cerr << endl
         << "Could not open ../OUTPUT/RS_mutation.dat"
         << endl; // Errore se il file non può essere aperto
    exit(3);
  }

  // Cicla sulle generazioni
  for (unsigned int j = 0; j < _ngen; j++) {
    double mean_loss_counter = 0.; // Contatore per la media del loss
    double best_loss = -1.;        // Miglior loss trovato in questa generazione

    field<Journey> parents = _journey_pop; // Copia la popolazione attuale come "genitori"
    field<Journey> sons(_journey_pop.size()); // Crea un field per i "figli"

    // Cicla su ogni individuo della popolazione
    for (unsigned int i = 0; i < _npop; i++) {
      Journey appo_journey = parents.at(i); // Prende un individuo padre

      // Applica gli operatori di mutazione con una data probabilità
      // SWAP MUTATION OPERATOR
      if (_rnd.Rannyu() < prob) {
        appo_journey.mutation_operator_swap();
      }

      // SWAP BLOCKS MUTATION OPERATOR
      if (_rnd.Rannyu() < prob) {
        appo_journey.mutation_operator_swap_blocks();
      }

      // SHIFT MUTATION OPERATOR
      if (_rnd.Rannyu() < prob) {
        appo_journey.mutation_operator_shift();
      }

      // FLIP MUTATION OPERATOR
      if (_rnd.Rannyu() < prob) {
        appo_journey.mutation_operator_flip_block();
      }

      appo_journey.compute_loss(_city_coords); // Ricalcola il loss per il figlio mutato
      sons.at(i) = appo_journey;               // Aggiunge il figlio al field dei figli
    }
    sons = sort_by_loss(sons); // Ordina i figli in base al loro loss

    // Selezione per la prossima generazione: metà dai figli migliori, metà dai genitori migliori
    uvec idx_sons = selection_operator((unsigned int)((_npop / 2) + _npop % 2));
    uvec idx_parents = selection_operator((unsigned int)(_npop / 2));

    // Riempie la nuova popolazione con i figli selezionati
    for (unsigned int k = 0; k < (int(_npop / 2) + _npop % 2); k++) {
      _journey_pop.at(k) = sons.at(idx_sons(k));
      _journey_pop.at(k).compute_loss(_city_coords); // Ricalcola il loss (anche se già fatto)
      double appo_loss = _journey_pop.at(k).getLoss();
      mean_loss_counter += appo_loss; // Aggiorna il contatore della media
      if (best_loss > appo_loss ||
          best_loss < 0.) { // Aggiorna il miglior loss della generazione
        best_loss = appo_loss;
      }
    }
    // Riempie la parte rimanente della nuova popolazione con i genitori selezionati
    for (unsigned int k = 0; k < (int(_npop / 2)); k++) {
      _journey_pop.at(k + int(_npop / 2) + _npop % 2) = parents.at(idx_parents(k));
      _journey_pop.at(k + int(_npop / 2) + _npop % 2)
          .compute_loss(_city_coords); // Ricalcola il loss
      double appo_loss = _journey_pop.at(k + int(_npop / 2) + _npop % 2).getLoss();
      mean_loss_counter += appo_loss; // Aggiorna il contatore della media
      if (best_loss > appo_loss || best_loss < 0.) {
        best_loss = appo_loss;
      }
    }
    double mean_loss =
        mean_loss_counter / _npop; // Calcola la media del loss per la generazione
    coutl << j << setw(12) << best_loss << setw(12)
          << mean_loss << endl; // Scrive i risultati nel file di output
  }
  coutl.close(); // Chiude il file
  _journey_pop = sort_by_loss(_journey_pop); // Riordina la popolazione finale
  return;
}

//////////////////////////////////////CROSSOVER///////////////////////////////////////////////////

// Metodo per applicare l'operatore di crossover a due genitori
void Population ::crossover(uvec idx_parents) {

  // Ottiene l'ordine delle città dai due genitori selezionati
  uvec parent1 = _journey_pop.at(idx_parents.at(0)).getCityOrder();
  uvec parent2 = _journey_pop.at(idx_parents.at(1)).getCityOrder();

  unsigned int M = parent1.size() - 1; // Lunghezza della parte "mobile" del percorso (esclusa la prima città)

  // Sceglie un punto di taglio casuale per il crossover (escludendo la prima città)
  unsigned int split_idx = static_cast<unsigned int>(_rnd.Rannyu(1, _ncities));

  // Estrae le parti "mancanti" (le code dopo il punto di taglio) dai genitori
  uvec missing_idx1 = parent1.subvec(split_idx, M);
  uvec missing_idx2 = parent2.subvec(split_idx, M);

  // Rimuove le code dai genitori per preparare il crossover
  parent1.shed_rows(split_idx, M);
  parent2.shed_rows(split_idx, M);

  // Creazione di maschere per trovare gli elementi delle "code" nell'altro genitore
  uvec mask1 = zeros<uvec>(
      _journey_pop.at(idx_parents.at(1)).getCityOrder().n_elem);
  uvec mask2 = zeros<uvec>(
      _journey_pop.at(idx_parents.at(0)).getCityOrder().n_elem);

  for (size_t i = 0;
       i < _journey_pop.at(idx_parents.at(1)).getCityOrder().n_elem; i++) {
    if (any(missing_idx1 == _journey_pop.at(idx_parents.at(1)).getCityOrder().at(
                                 i))) { // Se l'elemento di parent2 è presente in missing_idx1
      mask1(i) = 1; // Marca la posizione
    }
    if (any(missing_idx2 == _journey_pop.at(idx_parents.at(0)).getCityOrder().at(
                                 i))) { // Se l'elemento di parent1 è presente in missing_idx2
      mask2(i) = 1; // Marca la posizione
    }
  }
  // Trova gli indici degli elementi che corrispondono alla maschera
  uvec idx1 = find(mask1);
  uvec idx2 = find(mask2);

  // Crea il primo figlio unendo la testa del genitore1 con gli elementi di parent2 marcati da mask1
  uvec result1 = _journey_pop.at(idx_parents.at(1)).getCityOrder()(idx1);
  uvec son1 = join_vert(parent1, result1);

  // Crea il secondo figlio unendo la testa del genitore2 con gli elementi di parent1 marcati da mask2
  uvec result2 = _journey_pop.at(idx_parents.at(0)).getCityOrder()(idx2);
  uvec son2 = join_vert(parent2, result2);

  // Aggiorna l'ordine delle città dei genitori con i nuovi figli (i figli "sostituiscono" i genitori)
  _journey_pop.at(idx_parents.at(1)).setCityOrder(son1);
  _journey_pop.at(idx_parents.at(0)).setCityOrder(son2);
}

// Metodo per l'evoluzione della popolazione attraverso generazioni
void Population ::evolution(double p_c, double p_m1, double p_m2, double p_m3,
                           double p_m4) {
  // Apre il file di output specifico per l'evoluzione del continente corrente
  ofstream coutl("../OUTPUT/evolution" + to_string(_rank) + ".dat");
  if (coutl.is_open()) {
    coutl << "GEN:       BEST_LOSS:       AVERAGE_LOSS: "
          << endl; // Intestazione del file
  } else {
    cerr << endl
         << "Could not open ../OUTPUT/evolution" + to_string(_rank) + ".dat"
         << endl; // Errore se il file non può essere aperto
    exit(3);
  }
  _journey_pop = sort_by_loss(_journey_pop); // Ordina la popolazione iniziale

  // Ciclo sulle generazioni
  for (unsigned int j = 0; j < _ngen; j++) {
    double mean_loss_counter = 0.; // Contatore per la media del loss
    double best_loss = -1.;        // Miglior loss in questa generazione

    field<Journey> parents = _journey_pop; // Copia la popolazione attuale (genitori)

    // Cicla per creare nuove coppie di individui (metà della popolazione totale)
    for (unsigned int i = 0; i < _npop / 2; i++) {
      uvec parents_idx = selection_operator(2); // Seleziona due genitori
      if (_rnd.Rannyu() < p_c) { // Se la probabilità di crossover è soddisfatta
        crossover(parents_idx);  // Applica l'operatore di crossover
      } // Altrimenti i figli sono cloni dei genitori

      // Applicazione degli operatori di mutazione ai figli (o cloni)
      // SWAP MUTATION OPERATOR
      if (_rnd.Rannyu() < p_m1) {
        _journey_pop.at(parents_idx(0)).mutation_operator_swap();
      }
      if (_rnd.Rannyu() < p_m1) {
        _journey_pop.at(parents_idx(1)).mutation_operator_swap();
      }

      // SWAP BLOCKS MUTATION OPERATOR
      if (_rnd.Rannyu() < p_m2) {
        _journey_pop.at(parents_idx(0)).mutation_operator_swap_blocks();
      }
      if (_rnd.Rannyu() < p_m2) {
        _journey_pop.at(parents_idx(1)).mutation_operator_swap_blocks();
      }

      // SHIFT MUTATION OPERATOR
      if (_rnd.Rannyu() < p_m3) {
        _journey_pop.at(parents_idx(0)).mutation_operator_shift();
      }
      if (_rnd.Rannyu() < p_m3) {
        _journey_pop.at(parents_idx(1)).mutation_operator_shift();
      }

      // FLIP MUTATION OPERATOR
      if (_rnd.Rannyu() < p_m4) {
        _journey_pop.at(parents_idx(0)).mutation_operator_flip_block();
      }
      if (_rnd.Rannyu() < p_m4) {
        _journey_pop.at(parents_idx(1)).mutation_operator_flip_block();
      }

      // Ricalcola il loss per i due figli dopo mutazioni
      _journey_pop.at(parents_idx(1)).compute_loss(_city_coords);
      _journey_pop.at(parents_idx(0)).compute_loss(_city_coords);
    }
    _journey_pop = sort_by_loss(_journey_pop); // Riordina la popolazione aggiornata

    // Sostituisce una piccola percentuale della popolazione (i peggiori) con i migliori della generazione precedente (elitarismo)
    for (unsigned int k = 0; k < _journey_pop.size() / 5.; k++) {
      _journey_pop.at(k) = parents.at(k + 4 * static_cast<unsigned int>(_journey_pop.size() / 5.));
    }
    _journey_pop = sort_by_loss(_journey_pop); // Riordina nuovamente dopo l'elitarismo

    // Calcola il miglior loss e la media del loss per la generazione corrente
    for (unsigned int k = 0; k < _journey_pop.size(); k++) {
      double appo_loss = _journey_pop.at(k).getLoss();
      mean_loss_counter += appo_loss;
      if (appo_loss < best_loss || best_loss < 0.) { // Se è un nuovo miglior loss o il primo
        best_loss = appo_loss;
      }
    }
    double mean_loss = mean_loss_counter / double(_journey_pop.size());
    coutl << j << setw(12) << best_loss << setw(12)
          << mean_loss << endl; // Scrive i risultati nel file di output

    // Gestisce la migrazione tra i continenti a intervalli regolari
    if (j % _nmigr == 0 && j != 0 && _nContinents > 1) {
      for (unsigned int k = 0; k < _nContinents; k++) { // In media un paio di migrazioni per continente
        this->migration(); // Esegue la migrazione
      }
    }
  }
  coutl.close(); // Chiude il file di output
  return;
}

// Metodo per determinare il miglior percorso tra tutti i continenti (utilizzando MPI_Reduce)
void Population ::Olympics() {
  _journey_pop = sort_by_loss(_journey_pop); // Ordina la popolazione locale

  if (_nContinents > 1) { // Se ci sono più continenti (MPI)
    struct {
      double val; // Valore del loss
      int rank;   // Rank del processo
    } in, out;

    // `in.val` prende il loss dell'elemento con il loss maggiore nella popolazione locale (il peggiore se ordinato crescente)
    // Se la popolazione è ordinata in modo che il miglior loss sia all'ultimo indice, allora è corretto.
    // L'originale dice "loss massimo (o ultimo elemento, se già ordinato crescente)", il sort_by_loss ordina in "descend",
    // quindi l'elemento 0 ha il loss più alto e l'ultimo ha il loss più basso.
    // Se _journey_pop.at(_journey_pop.size()-1) è l'elemento con il loss più basso (il migliore), allora MPI_MINLOC è corretto.
    in.val = _journey_pop.at(_journey_pop.size() - 1).getLoss();
    in.rank = _rank;
    // Riduce tutti i valori di 'in' e trova il minimo valore ('MPI_MINLOC')
    // e il rank associato. Il risultato viene posto in 'out' sul processo con rank 0.
    MPI_Reduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
    if (_rank == 0) { // Solo il processo con rank 0 stampa l'ordine del vincitore
      print_order();
    }
  } else { // Se c'è solo un continente (non MPI)
    print_order(); // Stampa direttamente l'ordine del migliore localmente
  }
}

// Metodo per stampare l'ordine delle città del miglior percorso
void Population ::print_order() {
  _journey_pop = sort_by_loss(_journey_pop); // Assicura che la popolazione sia ordinata
  // Il miglior percorso è l'ultimo elemento dopo l'ordinamento "descend" (dal loss più grande al più piccolo)
  Journey best_journey = _journey_pop.at(_journey_pop.size() - 1);
  best_journey.compute_loss(_city_coords); // Ricalcola il loss del miglior percorso (sicurezza)
  _best_loss = best_journey.getLoss();     // Aggiorna il miglior loss globale

  // Se l'indice di simulazione è uguale al numero totale di passi (fine della simulazione)
  if (_sim_idx == _nblocks * _nsteps) {
    ofstream coutc("../OUTPUT/cities_order.dat"); // Apre il file per l'ordine delle città
    if (coutc.is_open()) {
      coutc << "CITIES ORDER:      ORDERED COORDINATES ---> X:     Y:"
            << endl; // Intestazione
      // Scrive l'ordine delle città e le loro coordinate
      for (unsigned int j = 0; j < best_journey.getCityOrder().size(); j++) {
        coutc << best_journey.getCityOrder().at(j)
              << setw(12)
              << _city_coords.at(best_journey.getCityOrder().at(j) - 1, 0)
              << setw(12)
              << _city_coords.at(best_journey.getCityOrder().at(j) - 1, 1)
              << endl;
      }
    } else {
      cerr << endl
           << "Couldn't open cities_order.dat, shutting down..."
           << endl; // Errore se il file non può essere aperto
      exit(4);
    }
    coutc.close(); // Chiude il file
  }
}

// Metodo per eseguire la simulazione utilizzando il metodo delle medie a blocchi
void Population ::blocks_mean() {

  Population POP; // Crea un nuovo oggetto Population per la simulazione
  POP.initialize(_nContinents,
                 _rank); // Inizializza il nuovo oggetto con i parametri

  double sum_average = 0, sum_ave2 = 0; // Variabili per il calcolo della media a blocchi

  // Ciclo sui blocchi
  for (unsigned int i = 0; i < _nblocks; i++) {
    double block_av = 0; // Media del blocco corrente

    // Ciclo sui passi all'interno di ogni blocco
    for (unsigned int j = 0; j < _nsteps; j++) {
      POP.create_pop(); // Crea una nuova popolazione per ogni step
      // Evolve la popolazione con i parametri di crossover e mutazione
      POP.evolution(0.7, 0.07, 0.07, 0.07, 0.07);
      POP._sim_idx++; // Incrementa l'indice della simulazione
      POP.Olympics(); // Esegue le "Olimpiadi" per trovare il miglior loss globale
      block_av += POP.getBestLoss(); // Aggiunge il miglior loss al contatore del blocco
    }
    if (_rank == 0) { // Solo il processo master scrive i risultati a blocchi
      ofstream coutf;
      double average;
      average = block_av / double(_nsteps); // Calcola la media del loss per il blocco
      sum_average += average;                // Aggiorna la somma delle medie dei blocchi
      sum_ave2 += average * average;         // Aggiorna la somma dei quadrati delle medie dei blocchi

      coutf.open("../OUTPUT/mb_best_loss.dat", ios::app); // Apre il file in modalità append
      coutf << setw(12) << i + 1 << setw(12) << average << setw(12)
            << sum_average / double(i + 1)
            << setw(12)
            << this->error(sum_average, sum_ave2, i + 1)
            << endl; // Scrive i risultati
      coutf.close(); // Chiude il file
    }
  }
  POP.finalize(); // Finalizza la simulazione
  return;
}

// Metodo per calcolare l'errore statistico (deviazione standard della media)
double Population ::error(double acc, double acc2, int blk) {
  if (blk <= 1)
    return 0.0; // Se c'è solo un blocco o meno, l'errore è 0
  else
    // Formula per la deviazione standard della media a blocchi
    return sqrt(fabs(acc2 / double(blk) - pow(acc / double(blk), 2)) /
                double(blk));
}

// Metodo per ottenere il miglior loss (la lunghezza del percorso più corta trovata)
double Population ::getBestLoss() { return _best_loss; }


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
