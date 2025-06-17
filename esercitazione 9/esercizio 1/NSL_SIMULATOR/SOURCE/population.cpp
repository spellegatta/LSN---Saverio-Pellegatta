/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cmath>    // Per funzioni matematiche come pow
#include <cstdlib>  // Per funzioni generali come exit
#include <string>   // Per l'uso delle stringhe
#include "population.h" // Include la definizione della classe Population
#include <cassert>  // Per la funzione assert

using namespace std;    // Utilizza lo spazio dei nomi standard
using namespace arma;   // Utilizza lo spazio dei nomi Armadillo

// Inizializza l'oggetto System leggendo i parametri dai file di input
void Population :: initialize(){

  int p1, p2; // Variabili per i numeri primi usati nell'inizializzazione del RNG
  // Apre il file Primes per leggere i numeri primi
  ifstream Primes("../INPUT/Primes");
  Primes >> p1 >> p2 ; // Legge i due numeri primi
  Primes.close(); // Chiude il file

  int seed[4]; // Array per il seed del generatore di numeri casuali
  // Apre il file seed.in per leggere il seed
  ifstream Seed("../INPUT/seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3]; // Legge i 4 valori del seed
  _rnd.SetRandom(seed,p1,p2); // Imposta il generatore di numeri casuali

  ifstream input("../INPUT/input.dat"); // Apre il file input.dat per leggere le configurazioni
  ofstream coutf; // Crea un oggetto ofstream per scrivere su file
  coutf.open("../OUTPUT/output.dat"); // Apre il file output.dat in modalità scrittura
  string property; // Stringa per leggere il nome della proprietà
  double delta; // Variabile non utilizzata, potrebbe essere un residuo
  // Cicla finché non si raggiunge la fine del file di input
  while ( !input.eof() ){
    input >> property; // Legge il nome della proprietà
    if( property == "NCITIES" ){
      input >> _ncities; // Legge il numero di città
    }else if( property == "SIDE_RADIUS" ){
      input >> _side_radius; // Legge il lato del quadrato o il raggio del cerchio
    }else if( property == "SIMULATION_TYPE" ){
      input >> _sim_type; // Legge il tipo di simulazione (0 per cerchio, 1 per quadrato)
      if(_sim_type > 1){
        cerr << "PROBLEM: unknown simulation type" << endl; // Errore per tipo di simulazione sconosciuto
        exit(EXIT_FAILURE); // Esce con errore
      }
      if(_sim_type == 0){
        coutf << "CITIES ON A CIRCULAR CONFIGURATION"  << endl; // Stampa il tipo di configurazione
        _city_coords.set_size(_ncities, 2); // Ridimensiona la matrice delle coordinate
        ofstream coutc; // Crea un oggetto ofstream
        coutc.open("../OUTPUT/coordinates.dat"); // Apre il file coordinates.dat
        coutc << "X:     Y: " << endl; // Scrive l'intestazione
        for (unsigned int j=0; j<_ncities; j++){ // Genera le coordinate delle città su un cerchio
          double theta=_rnd.Rannyu(0, 2*M_PI); // Genera un angolo casuale
          _city_coords(j,1)=_side_radius*sin(theta); // Calcola la coordinata Y
          _city_coords(j, 0)=_side_radius*cos(theta); // Calcola la coordinata X
          coutc << _city_coords(j,0) << setw(12) << _city_coords(j, 1) << endl; // Scrive le coordinate
        }
        coutc.close(); // Chiude il file delle coordinate
      }
      else if(_sim_type == 1) {
        _city_coords.set_size(_ncities, 2); // Ridimensiona la matrice delle coordinate
        coutf << "CITIES IN A SQUARE AREA" << endl; // Stampa il tipo di configurazione
        ofstream coutc; // Crea un oggetto ofstream
        coutc.open("../OUTPUT/coordinates.dat"); // Apre il file coordinates.dat
        coutc << "X:     Y: " << endl; // Scrive l'intestazione
        for (unsigned int j=0; j<_ncities; j++){ // Genera le coordinate delle città in un quadrato
          _city_coords(j,1)=_rnd.Rannyu(0, _side_radius); // Calcola la coordinata Y
          _city_coords(j, 0)=_rnd.Rannyu(0, _side_radius); // Calcola la coordinata X
          coutc << _city_coords(j,1) << setw(12) << _city_coords(j, 0) << endl; // Scrive le coordinate
        }
      }
    } else if( property == "NPOP" ){
      input >> _npop; // Legge la dimensione della popolazione
    } else if( property == "NGEN" ){
      input >> _ngen; // Legge il numero di generazioni
    }else if( property == "SELECTION_P" ){
      input >> _p; // Legge l'esponente per l'operatore di selezione
    }else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << endl; // Messaggio di completamento lettura input
      break; // Esce dal ciclo
    } else cerr << "PROBLEM: unknown input" << endl; // Errore per input sconosciuto
  }
  input.close(); // Chiude il file di input
  _journey_pop.set_size(_npop); // Ridimensiona il campo della popolazione
  create_pop(); // Crea la popolazione iniziale
  _journey_pop = sort_by_loss(_journey_pop); // Ordina la popolazione in base alla "perdita"
  coutf << "System initialized!" << endl; // Messaggio di inizializzazione completata
  coutf.close(); // Chiude il file di output
  return;
}

// Crea la popolazione iniziale di percorsi casuali
void Population :: create_pop(){
  for (unsigned int i=0; i<_npop; i++){ // Itera per il numero di individui nella popolazione
    Journey appo_journey(_ncities, _rnd, _city_coords); // Crea un nuovo percorso casuale
    _journey_pop.at(i)=appo_journey; // Aggiunge il percorso alla popolazione
  }
  return;
}

// Ordina una popolazione di percorsi in base alla loro "perdita" (dal peggiore al migliore)
field<Journey> Population :: sort_by_loss(field<Journey> field_pop){
  vec appo_losses(field_pop.size()); // Vettore temporaneo per memorizzare le perdite
  for(unsigned int i=0; i<field_pop.size(); i++){ // Itera su ogni percorso nella popolazione
    appo_losses.at(i)=field_pop.at(i).getLoss(); // Ottiene la perdita del percorso
  }
  uvec idx = sort_index(appo_losses, "descend"); // Ottiene gli indici di ordinamento in ordine decrescente
  field<Journey> appo(_npop); // Crea un nuovo campo per la popolazione ordinata
  for(unsigned int j=0; j<field_pop.size(); j++){ // Riordina i percorsi in base agli indici
    appo.at(j)=field_pop.at(idx.at(j));
  }
  field_pop=appo; // Assegna la popolazione ordinata
  return field_pop; // Restituisce la popolazione ordinata
}

// Finalizza la simulazione
void Population :: finalize(){
  _rnd.SaveSeed(); // Salva il seed del generatore di numeri casuali
  ofstream coutf; // Crea un oggetto ofstream
  coutf.open("../OUTPUT/output.dat",ios::app); // Apre il file output.dat in modalità append
  coutf << "Simulation completed!" << endl; // Scrive il messaggio di completamento
  coutf.close(); // Chiude il file
  return;
}

/////////////////////////////OPERATORE DI SELEZIONE///////////////////////////////////////////////////

// Operatore di selezione: seleziona gli individui per la riproduzione
uvec Population :: selection_operator(unsigned int n_people){
  uvec idx(n_people); // Vettore per memorizzare gli indici degli individui selezionati
  for(unsigned int i=0; i<n_people; i++){
    // Seleziona un indice in base a una distribuzione di potenza (_p)
    // Questo favorisce la selezione di individui con "perdita" inferiore (migliori)
    idx.at(i)=(unsigned int)(_npop*pow(_rnd.Rannyu(), _p));
  }
  return idx; // Restituisce gli indici selezionati
}

/////////////////////////////////////RICERCA CASUALE//////////////////////////////////////////////////
// Esegue una ricerca casuale applicando solo operatori di mutazione
void Population :: rand_search_mutation(double prob){
  ofstream coutl("../OUTPUT/RS_mutation.dat"); // Apre il file per i risultati della mutazione
  if(coutl.is_open()){
    coutl << "GEN:       BEST_LOSS:       AVERAGE_LOSS: " << endl; // Scrive l'intestazione
  } else{
    cerr << endl << "Could not open ../OUTPUT/RS_mutation.dat" << endl; // Errore se il file non si apre
    exit(3); // Esce con errore
  }
  for(unsigned int j=0; j<_ngen; j++){ // Cicla sulle generazioni
    double mean_loss_counter=0.; // Contatore per la perdita media
    double best_loss=-1.; // Migliore perdita trovata
    field<Journey> parents=_journey_pop; // Copia la popolazione attuale (genitori)
    field<Journey> sons(_journey_pop.size()); // Crea un campo per i figli

    for(unsigned int i=0; i<_npop; i++){ // Cicla sugli individui della popolazione
      Journey appo_journey=parents.at(i); // Copia il percorso del genitore

      // SWAP MUTATION OPERATOR
      if(_rnd.Rannyu()<prob){ // Applica la mutazione con una data probabilità
        appo_journey.mutation_operator_swap();
      }

      // SWAP BLOCKS MUTATION OPERATOR
      if(_rnd.Rannyu()<prob){ // Applica la mutazione con una data probabilità
        appo_journey.mutation_operator_swap_blocks();
      }

      // SHIFT MUTATION OPERATOR
      if(_rnd.Rannyu()<prob){ // Applica la mutazione con una data probabilità
        appo_journey.mutation_operator_shift();
      }

      // FLIP MUTATION OPERATOR
      if(_rnd.Rannyu()<prob){ // Applica la mutazione con una data probabilità
        appo_journey.mutation_operator_flip_block();
      }

      appo_journey.compute_loss(_city_coords); // Ricalcola la perdita per il figlio mutato
      sons.at(i)=appo_journey; // Aggiunge il figlio alla popolazione dei figli
    }
    sons=sort_by_loss(sons); // Ordina i figli per perdita

    // Selezione di figli e genitori per formare la nuova generazione
    uvec idx_sons = selection_operator((unsigned int)((_npop/2)+_npop%2)); // Seleziona figli
    uvec idx_parents = selection_operator((unsigned int)(_npop/2)); // Seleziona genitori
    for(unsigned int k=0; k<(int(_npop/2)+_npop%2); k++){
      _journey_pop.at(k)=sons.at(idx_sons(k)); // Sostituisce i percorsi con i figli selezionati
      _journey_pop.at(k).compute_loss(_city_coords); // Ricalcola la perdita
      double appo_loss=_journey_pop.at(k).getLoss(); // Ottiene la perdita
      mean_loss_counter+=appo_loss; // Aggiorna il contatore della perdita media
      if(best_loss>appo_loss or best_loss<0.){ // Aggiorna la migliore perdita
        best_loss=appo_loss;
      }
    }
    for(unsigned int k=0; k<(int(_npop/2)); k++){
      _journey_pop.at(k+int(_npop/2)+_npop%2)=parents.at(idx_parents(k)); // Sostituisce i percorsi con i genitori selezionati
      _journey_pop.at(k+int(_npop/2)+_npop%2).compute_loss(_city_coords); // Ricalcola la perdita
      double appo_loss=_journey_pop.at(k+int(_npop/2)+_npop%2).getLoss(); // Ottiene la perdita
      mean_loss_counter+=appo_loss; // Aggiorna il contatore della perdita media
      if(best_loss>appo_loss or best_loss<0.){ // Aggiorna la migliore perdita
        best_loss=appo_loss;
      }
    }
    double mean_loss=mean_loss_counter/_npop; // Calcola la perdita media della generazione
    coutl << j << setw(12) << best_loss << setw(12) << mean_loss << endl; // Scrive i risultati sul file
  }
  coutl.close(); // Chiude il file
  _journey_pop=sort_by_loss(_journey_pop); // Riordina la popolazione finale
  return;
}

//////////////////////////////////////CROSSOVER///////////////////////////////////////////////////

// Operatore di crossover: combina due percorsi per crearne due nuovi
void Population :: crossover(uvec idx_parents){

  uvec parent1=_journey_pop.at(idx_parents.at(0)).getCityOrder(); // Ottiene l'ordine delle città del primo genitore
  uvec parent2=_journey_pop.at(idx_parents.at(1)).getCityOrder(); // Ottiene l'ordine delle città del secondo genitore

  unsigned int M = parent1.size()-1; // Dimensione della parte mobile del percorso

  // Sceglie un punto di divisione casuale per il crossover
  unsigned int split_idx= static_cast<unsigned int>(_rnd.Rannyu(1, _ncities));

  // Estrae la parte finale dei percorsi che verrà scambiata
  uvec missing_idx1=parent1.subvec(split_idx, M);
  uvec missing_idx2=parent2.subvec(split_idx, M);

  // Rimuove la parte finale dai genitori
  parent1.shed_rows(split_idx, M);
  parent2.shed_rows(split_idx, M);
  
  // Crea maschere per identificare le città mancanti nei nuovi figli
  uvec mask1 = zeros<uvec>(_journey_pop.at(idx_parents.at(1)).getCityOrder().n_elem);
  uvec mask2 = zeros<uvec>(_journey_pop.at(idx_parents.at(0)).getCityOrder().n_elem);

  for (size_t i = 0; i < _journey_pop.at(idx_parents.at(1)).getCityOrder().n_elem; i++) {
    // Se la città è presente nella parte mancante dell'altro genitore, segnala nella maschera
    if (any(missing_idx1 == _journey_pop.at(idx_parents.at(1)).getCityOrder().at(i))) {
      mask1(i) = 1;
    }
    if (any(missing_idx2 == _journey_pop.at(idx_parents.at(0)).getCityOrder().at(i))) {
      mask2(i) = 1;
    }
  }
  // Trova gli indici delle città da inserire per completare i figli
  uvec idx1  = find(mask1);
  uvec idx2  = find(mask2);

  // Costruisce il primo figlio unendo la parte iniziale del genitore 1 con le città ordinate del genitore 2
  uvec result1 = _journey_pop.at(idx_parents.at(1)).getCityOrder()(idx1);
  uvec son1=join_vert(parent1, result1);

  // Costruisce il secondo figlio unendo la parte iniziale del genitore 2 con le città ordinate del genitore 1
  uvec result2 = _journey_pop.at(idx_parents.at(0)).getCityOrder()(idx2);
  uvec son2=join_vert(parent2, result2);

  // Aggiorna i percorsi dei genitori con i nuovi figli generati
  _journey_pop.at(idx_parents.at(1)).setCityOrder(son1);
  _journey_pop.at(idx_parents.at(0)).setCityOrder(son2);

}

// Esegue il processo di evoluzione della popolazione (algoritmo genetico)
void Population :: evolution(double p_c, double p_m1, double p_m2, double p_m3, double p_m4){
  ofstream coutl("../OUTPUT/evolution.dat"); // Apre il file per i risultati dell'evoluzione
  if(coutl.is_open()){
    coutl << "GEN:       BEST_LOSS:       AVERAGE_LOSS: " << endl; // Scrive l'intestazione
  } else{
    cerr << endl << "Could not open ../OUTPUT/evolution.dat" << endl; // Errore se il file non si apre
    exit(3); // Esce con errore
  }
  _journey_pop=sort_by_loss(_journey_pop); // Ordina la popolazione iniziale

  for(unsigned int j=0; j<_ngen; j++){ // Cicla sulle generazioni
    double mean_loss_counter=0.; // Contatore per la perdita media
    double best_loss=-1.; // Migliore perdita trovata
    field<Journey> parents=_journey_pop; // Copia la popolazione attuale (genitori)

    for(unsigned int i=0; i<_npop/2; i++){ // Cicla per metà della popolazione (generiamo coppie di figli)
      uvec parents_idx=selection_operator(2); // Seleziona due genitori
      if(_rnd.Rannyu()<p_c){ // Se una condizione casuale è soddisfatta (probabilità di crossover)
        crossover(parents_idx); // Applica l'operatore di crossover
      } // Se no, i figli sono cloni dei genitori

      // Applica vari operatori di mutazione con le rispettive probabilità
      // SWAP MUTATION OPERATOR
      if(_rnd.Rannyu()<p_m1){
        _journey_pop.at(parents_idx(0)).mutation_operator_swap();
      }
      if(_rnd.Rannyu()<p_m1){
        _journey_pop.at(parents_idx(1)).mutation_operator_swap();
      }

      // SWAP BLOCKS MUTATION OPERATOR
      if(_rnd.Rannyu()<p_m2){
        _journey_pop.at(parents_idx(0)).mutation_operator_swap_blocks();
      }
      if(_rnd.Rannyu()<p_m2){
        _journey_pop.at(parents_idx(1)).mutation_operator_swap_blocks();
      }

      // SHIFT MUTATION OPERATOR
      if(_rnd.Rannyu()<p_m3){
        _journey_pop.at(parents_idx(0)).mutation_operator_shift();
      }
      if(_rnd.Rannyu()<p_m3){
        _journey_pop.at(parents_idx(1)).mutation_operator_shift();
      }

      // FLIP MUTATION OPERATOR
      if(_rnd.Rannyu()<p_m4){
        _journey_pop.at(parents_idx(0)).mutation_operator_flip_block();
      }
      if(_rnd.Rannyu()<p_m4){
        _journey_pop.at(parents_idx(1)).mutation_operator_flip_block();
      }

      // Ricalcola la perdita per i due figli/individui mutati
      _journey_pop.at(parents_idx(1)).compute_loss(_city_coords);
      _journey_pop.at(parents_idx(0)).compute_loss(_city_coords);

    }
    _journey_pop=sort_by_loss(_journey_pop); // Riordina la popolazione dopo mutazioni e crossover
    // Sostituisce la metà peggiore della popolazione con la metà migliore dei genitori originali
    for(unsigned int k=0; k<_journey_pop.size()/2; k++){
      _journey_pop.at(k)=parents.at(k+static_cast<unsigned int>(_journey_pop.size()/2));
    }
    _journey_pop=sort_by_loss(_journey_pop); // Riordina nuovamente la popolazione aggiornata

    // Calcola la migliore perdita e la perdita media per la generazione corrente
    for (unsigned int k=0; k<_journey_pop.size(); k++){
      double appo_loss=_journey_pop.at(k).getLoss();
      mean_loss_counter+=appo_loss;
      if(appo_loss<best_loss or best_loss<0.){
        best_loss=appo_loss;
      }
    }
    double mean_loss=mean_loss_counter/double(_journey_pop.size()); // Calcola la perdita media
    coutl << j << setw(12) << best_loss << setw(12) << mean_loss << endl; // Scrive i risultati sul file
  }
  coutl.close(); // Chiude il file
  return;
}

// Stampa l'ordine delle città del miglior percorso trovato
void Population :: print_order(){
  _journey_pop=sort_by_loss(_journey_pop); // Assicura che la popolazione sia ordinata
  Journey best_journey= _journey_pop.at(_journey_pop.size()-1); // Ottiene il miglior percorso (ultima posizione dopo l'ordinamento "descend")
  ofstream coutc("../OUTPUT/cities_order.dat"); // Apre il file per scrivere l'ordine delle città
  if (coutc.is_open()){
    coutc << "CITIES ORDER:       ORDERED COORDINATES ---> X:      Y:" << endl; // Scrive l'intestazione
    for(unsigned int j=0; j< best_journey.getCityOrder().size(); j++){
      // Scrive il numero della città e le sue coordinate ordinate
      coutc << best_journey.getCityOrder().at(j) <<setw(12) << _city_coords.at(best_journey.getCityOrder().at(j)-1, 0) << setw(12) << _city_coords.at(best_journey.getCityOrder().at(j)-1, 1) <<  endl;
    }
  } else{
    cerr << endl << "Couldn't open cities_order.dat, shutting down..." << endl; // Errore se il file non si apre
    exit(4); // Esce con errore
  }
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
