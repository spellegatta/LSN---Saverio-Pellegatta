/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream> // Per input e output da console
#include <math.h>   // Per funzioni matematiche come sqrt, sin, cos
#include "journey.h"// Include la definizione della classe Journey
#include <cassert>  // Per la funzione assert, usata per debugging e validazione

using namespace std; // Utilizza lo spazio dei nomi standard

// Costruttore della classe Journey
Journey :: Journey(unsigned int ncities, Random rnd, mat coord){
    // Inizializza il numero di città con il valore passato
    _ncities=ncities;
    // Assegna il generatore di numeri casuali
    _rnd=rnd;
    // Ridimensiona il vettore che conterrà l'ordine delle città
    _cities_order.set_size(_ncities);
    // Crea un vettore base con i numeri delle città da 1 a ncities
    uvec base = regspace<uvec>(1, _ncities);
    // Estrae una sottosezione del vettore base, escludendo la prima città (che rimane fissa)
    uvec appo=base.subvec(1, base.size()-1);
    // Mescola casualmente le città rimanenti
    appo=shuffle(appo);
    // Imposta la prima città nell'ordine come la città di partenza (assunta essere la città 1)
    _cities_order.at(0)=base.at(0);
    // Assegna le città mescolate alla parte rimanente del percorso
    _cities_order.subvec(1, _cities_order.size()-1)=appo;
    // Controlla la validità del percorso appena creato
    check_validity();
    // Calcola la "perdita" (o costo/lunghezza) del percorso
    compute_loss(coord);
}

// Costruttore di default della classe Journey
Journey::Journey(){
    // Inizializza il numero di città a zero
    _ncities=0;
    // Crea un'istanza del generatore di numeri casuali
    Random rnd;
    _rnd=rnd;
    // Inizializza il vettore dell'ordine delle città con dimensione zero
    _cities_order.set_size(_ncities);
    // Inizializza la perdita a un valore negativo per indicare che non è stata calcolata
    _loss=-1.;
}

// Funzione per calcolare la lunghezza totale del percorso (loss)
void Journey :: compute_loss(mat coord){
    _loss=0; // Inizializza la perdita a zero
    // Itera su tutte le città del percorso
    for(unsigned int j=0; j< coord.n_rows; j++){
        // Calcola la differenza tra le coordinate della città corrente e della città successiva
        // Utilizza l'operatore modulo per collegare l'ultima città alla prima (percorso circolare)
        rowvec diff = coord.row(_cities_order.at(j)-1)-coord.row((_cities_order.at((j+1)%_cities_order.size()))-1);
        // Aggiunge la norma euclidea (distanza) tra le due città alla perdita totale
        _loss += norm(diff, 2);
    }
}

// Funzione per ottenere il valore della perdita del percorso
double Journey :: getLoss(){
    return _loss; // Restituisce la perdita
}

// Funzione per ottenere l'ordine delle città nel percorso
uvec Journey :: getCityOrder(){
    return _cities_order; // Restituisce il vettore con l'ordine delle città
}

// Funzione per impostare l'ordine delle città nel percorso
void Journey :: setCityOrder(uvec cities_order){
    _cities_order=cities_order; // Assegna il nuovo ordine delle città
}

// Funzione per verificare la validità del percorso
void Journey :: check_validity(){
    // Assicura che il viaggio inizi sempre dalla città 1
    assert(("The salesman should start his journey in city 1",_cities_order.at(0)==1));

    // Controlla che nessuna città sia stata visitata due volte
    for(unsigned int i=0; i<_cities_order.size(); i++){
        for(unsigned int j=i+1; j<_cities_order.size(); j++){ // Confronta ogni elemento solo una volta per evitare doppioni e auto-confronti
            if(_cities_order.at(i)==_cities_order.at(j)){
                cerr << endl << "A city has been visited twice, exiting..." << endl;
                exit(1); // Esce dal programma con codice di errore
            }
        }
    }

    // Controlla che tutte le città siano state visitate
    for(unsigned int i=0; i < _ncities; i++){
        // Crea una maschera booleana: 1 se il valore in _cities_order è uguale a (i+1), 0 altrimenti
        uvec mask= _cities_order==(i+1);
        // Se 'any(mask)' restituisce false (cioè la maschera è tutta zeri), significa che la città (i+1) non è stata visitata
        if (!any(mask)){
            cerr << endl << "At least a city has not been visited" << endl;
            exit(2); // Esce dal programma con codice di errore
        }
    }
    return; // Il percorso è valido
}

///////////////////////////OPERATORI DI MUTAZIONE////////////////////////////////////////////////

// Operatore di mutazione: scambia due città casuali nel percorso (escludendo la prima)
void Journey :: mutation_operator_swap(){
    unsigned int i, j;
    // Genera due indici casuali (escludendo l'indice 0, che corrisponde alla città di partenza fissa)
    i=int(_rnd.Rannyu(1, _ncities));
    j=int(_rnd.Rannyu(1, _ncities));

    // Assicura che gli indici generati siano positivi
    assert(("City order vector indexes should be positive", i>0 and j>0));

    // Scambia le città agli indici i e j
    unsigned int appo=_cities_order.at(i);
    _cities_order.at(i)=_cities_order.at(j);
    _cities_order.at(j)=appo;

    check_validity(); // Controlla la validità del percorso dopo la mutazione
    return;
}

// Operatore di mutazione: sposta un blocco di città
void Journey::mutation_operator_shift() {
    // 1) Parametri e setup
    const unsigned int N = _cities_order.size();           // Dimensione totale del percorso (inclusa la città fissa)
    const unsigned int M = N - 1;                           // Lunghezza della parte "mobile" del percorso (esclusa la città di partenza)
    uvec pos0            = _cities_order.subvec(0,0);       // Vettore contenente solo la prima città (fissa)
    uvec mobile          = _cities_order.subvec(1, N-1);    // Vettore contenente la parte modificabile del percorso

    // 2) Estrai parametri casuali
    unsigned int i     = unsigned(_rnd.Rannyu(0, M));           // Indice di inizio del blocco da spostare [0..M-1]
    unsigned int m     = 1 + unsigned(_rnd.Rannyu(0, M-1));     // Lunghezza del blocco da spostare [1..M]
    // Indice di destinazione per l'inserimento del blocco, con wrap-around
    unsigned int f     = (i + unsigned(_rnd.Rannyu(0, M))) % M;

    // 3) Prendi il blocco (anche se "wrap" oltre la fine del vettore)
    uvec block_vec(m); // Vettore temporaneo per memorizzare il blocco
    for (unsigned k = 0; k < m; ++k) {
        block_vec(k) = mobile[(i + k) % M]; // Prende gli elementi con wrap-around
    }

    // 4) Prepara il nuovo vettore "mobile" ricostruito
    uvec new_mobile(M); // Nuovo vettore mobile
    std::vector<bool> filled(M,false); // Vettore booleano per tenere traccia delle posizioni riempite

    // 4a) Inserisci il blocco alla destinazione, con wrap-around
    for (unsigned k = 0; k < m; ++k) {
        unsigned pos = (f + k) % M; // Calcola la posizione con wrap-around
        new_mobile(pos) = block_vec(k); // Inserisce l'elemento del blocco
        filled[pos]     = true;         // Segna la posizione come riempita
    }

    // 4b) Rimetti "in ordine" tutti gli altri elementi
    //     partendo subito dopo il blocco estratto nella vecchia mobile
    unsigned rest_idx = (i + m) % M; // Primo elemento "non-bloccato" dopo l'estrazione
    for (unsigned pos = 0; pos < M; ++pos) {
        if (!filled[pos]) { // Se la posizione non è ancora stata riempita
            new_mobile(pos) = mobile[rest_idx]; // Inserisce un elemento dal resto della vecchia mobile
            rest_idx = (rest_idx + 1) % M;      // Avanza all'elemento successivo nel resto
        }
    }

    // 5) Ricostruisci l'ordine completo
    _cities_order.set_size(N); // Ridimensiona il vettore dell'ordine delle città
    _cities_order[0] = pos0(0); // Reinserisce la prima città fissa
    _cities_order.subvec(1, N-1) = new_mobile; // Assegna il nuovo vettore mobile

    // 6) Check finale
    check_validity(); // Controlla la validità del percorso dopo la mutazione
}

// Operatore di mutazione: scambia due blocchi di città
void Journey::mutation_operator_swap_blocks() {
    // 1) dimensioni
    unsigned int N = static_cast<unsigned int>(_cities_order.n_elem); // Numero totale di elementi nel percorso
    unsigned int M = N - 1;                                            // Lunghezza della parte "mobile"

    // 2) view diretto sulla subvec [1..N-1]
    uvec mobile = _cities_order.subvec(1, N-1); // Ottiene la parte mobile del percorso

    // 3) scegli lunghezza del blocco in [1..M/2]
    unsigned int block_length = 1 + static_cast<unsigned int>(_rnd.Rannyu(0, M/2)); // Lunghezza del blocco da scambiare

    // 4) scegli start in [0..M-1]
    unsigned int start1 = static_cast<unsigned int>(_rnd.Rannyu(0, M)); // Indice di inizio del primo blocco

    // 5) scegli di quanto spostare in [1..M-1]
    unsigned int nshifts = 1 + static_cast<unsigned int>(_rnd.Rannyu(0, M-1)); // Numero di posizioni per il secondo blocco

    // 6) fai lo swap ciclico con modulo M
    for (unsigned int k = 0; k < block_length; ++k) {
        unsigned int idx1 = (start1 + k) % M;                  // Indice del primo elemento del blocco 1 (con wrap)
        unsigned int idx2 = (start1 + nshifts + k) % M;        // Indice del primo elemento del blocco 2 (con wrap)
        swap(mobile(idx1), mobile(idx2));                       // Scambia gli elementi
    }

    // 7) la subvec modifica direttamente _cities_order, quindi basta validare
    check_validity(); // Controlla la validità del percorso dopo la mutazione
}

// Operatore di mutazione: inverte l'ordine di un blocco di città
void Journey :: mutation_operator_flip_block(){
    // Genera una lunghezza casuale per il blocco da invertire (escludendo la prima città)
    unsigned int block_length=(unsigned int)(_rnd.Rannyu(1, int(_cities_order.size()-1)));
    // Ottiene la parte mobile del percorso
    uvec mobile=_cities_order.subvec(1, _cities_order.size()-1);
    // Ottiene la prima città (fissa)
    uvec pos0=_cities_order.subvec(0, 0);

    // Genera un indice di partenza casuale per il blocco all'interno della parte mobile
    unsigned int start = (unsigned int)(_rnd.Rannyu(0,mobile.size()));

    // Caso 1: il blocco non "wrappa" (non supera la fine del vettore)
    if(start+block_length<=mobile.size()){
        uvec appo=mobile.subvec(start, start+block_length-1); // Estrae il blocco
        appo=flipud(appo);                                     // Inverte l'ordine del blocco
        mobile.subvec(start, start+block_length-1)=appo;       // Reinserisce il blocco invertito
    } else{ // Caso 2: il blocco "wrappa" (supera la fine del vettore)
        // Estrae la prima parte del blocco (dall'inizio al wraparound)
        uvec appo1=mobile.subvec(start, mobile.size()-1);
        // Estrae la seconda parte del blocco (dopo il wraparound, dall'inizio del vettore)
        uvec appo2=mobile.subvec(0, start+block_length-mobile.size()-1);
        // Unisce le due parti e inverte l'intero blocco
        uvec appo=flipud(join_vert(appo1, appo2));
        // Reinserisce la prima parte del blocco invertito
        mobile.subvec(start, mobile.size()-1)=appo.subvec(0, mobile.size()-1-start);
        // Reinserisce la seconda parte del blocco invertito
        mobile.subvec(0, start+block_length-mobile.size()-1)=appo.subvec(mobile.size()-start, appo.size()-1);
    }
    // Ricostruisce l'ordine completo delle città unendo la città fissa con la parte mobile modificata
    _cities_order=join_vert(pos0, mobile);

    check_validity(); // Controlla la validità del percorso dopo la mutazione
}

// Operatore di assegnazione per la classe Journey
Journey& Journey::operator=(const Journey& other) {
    // Il commento seguente è stato lasciato nel codice originale, probabilmente per debugging.
    // cout << endl <<"before cities " << other._ncities << " order   " << other._cities_order.t() << "   loss " << other._loss << endl;

    // Controlla per auto-assegnazione per evitare problemi
    if (this != &other) {
        _ncities       = other._ncities;        // Copia il numero di città
        _rnd           = other._rnd;            // Copia il generatore di numeri casuali
        _cities_order = other._cities_order;   // Copia l'ordine delle città
        _loss          = other._loss;           // Copia la perdita del percorso
    }
    // Il commento seguente è stato lasciato nel codice originale, probabilmente per debugging.
    // cout << endl <<"after cities " << _ncities << " order   " << _cities_order.t() << "   loss " << _loss << endl;
    return *this; // Restituisce un riferimento all'oggetto corrente
}