// functions.cpp

#include <vector>
#include <iostream>
#include "random.h"
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

// Calcola la media di nThrows numeri estratti con distribuzione uniforme [0,1)
double UniformMean(Random* rnd, unsigned int nThrows){
    double cumsum = 0;                 // somma cumulativa dei numeri estratti
    for (unsigned int i = 0; i < nThrows; i++) {
        cumsum += rnd->Rannyu();       // estrazione uniforme
    }
    return cumsum / nThrows;           // media del campione
}

// Calcola la media di nThrows numeri estratti con distribuzione esponenziale
double ExponentialMean(Random* rnd, unsigned int nThrows, double lambda){
    double cumsum = 0;                 // somma cumulativa delle estrazioni
    for (unsigned int i = 0; i < nThrows; i++) {
        cumsum += rnd->Exponential(lambda);  // estrazione esponenziale (lambda)
    }
    return cumsum / nThrows;           // media del campione
}

// Calcola la media di nThrows numeri estratti con distribuzione di Cauchy-Lorentz
double CauchyLorentzMean(Random* rnd, unsigned int nThrows, double mu, double gamma){
    double cumsum = 0;                 // somma cumulativa delle estrazioni
    for (unsigned int i = 0; i < nThrows; i++) {
        cumsum += rnd->CauchyLorentz(mu, gamma);  // estrazione Cauchy-Lorentz (mu, gamma)
    }
    return cumsum / nThrows;           // media del campione
}

// Scrive su file tutti i valori del vettore vec, uno per riga
void outData(vector<double> vec, const char file[]){
    ofstream output;
    output.open(file);
    if (output.fail()) {
        cout << endl
             << "Errore di apertura dello stream di output, "
                "il file su cui vuoi trascrivere è nella stessa cartella dell'eseguibile?"
             << endl;
        return;                       // esco se il file non si apre
    }
    for (unsigned int i = 0; i < vec.size(); i++) {
        output << vec.at(i) << endl;  // stampo elemento i-esimo
    }
    output.close();                   // chiudo lo stream
}
