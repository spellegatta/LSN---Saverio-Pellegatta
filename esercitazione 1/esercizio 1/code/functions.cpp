// functions.cpp

#include <vector>
#include <iostream>
#include "random.h"
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

// Calcola la media di un singolo blocco di 'throws' estrazioni uniformi [0,1)
double OneBlockMean(Random* rnd, unsigned int throws){
    double x;                     // variabile per il valore estratto
    double cum_sum = 0;           // somma cumulativa dei valori estratti
    for (unsigned int i = 0; i < throws; i++) {
        x = rnd->Rannyu();        // estraggo un numero casuale
        cum_sum += x;             // accumulo nella somma
    }
    return cum_sum / throws;      // restituisco la media del blocco
}

// Scrive su file i vettori mean e StdDev in due colonne
void outData(std::vector<double> mean, std::vector<double> StdDev, const char file[]){
    ofstream output;
    output.open(file);
    if (output.fail()) {
        cout << endl
             << "Errore di apertura dello stream di output, "
                "il file su cui vuoi trascrivere è nella stessa cartella dell'eseguibile?"
             << endl;
        return;                   // esco se il file non si apre
    }
    for (unsigned int i = 0; i < mean.size(); i++) {
        output << mean.at(i)
               << setw(12)       // allineamento della seconda colonna
               << StdDev.at(i)
               << endl;
    }
    output.close();             // chiudo lo stream
}

// Costruisce progressivamente la media e la sua deviazione standard per blocchi crescenti
void MeanAndDevStd(Random* rnd, unsigned int nBlocksMax, unsigned int nThrows){
    double CumSum = 0;             // somma cumulativa delle medie dei blocchi
    double SquaredCumSum = 0;      // somma cumulativa dei quadrati delle medie
    double x;                      // variabile di appoggio
    vector<double> mean, devstd;   // vettori di output: media cumulativa e dev. std.

    // Caso del primo blocco: deviazione standard = 0 per definizione
    x = OneBlockMean(rnd, nThrows);
    CumSum += x;
    SquaredCumSum += pow(x, 2);
    mean.push_back(CumSum);        // media cumulativa con 1 solo blocco
    devstd.push_back(0);           // dev. std. nulla

    // Ciclo per i blocchi da 2 a nBlocksMax
    for (unsigned int i = 2; i <= nBlocksMax; i++) {
        x = OneBlockMean(rnd, nThrows);
        CumSum += x;
        SquaredCumSum += pow(x, 2);
        mean.push_back(CumSum / double(i));  
        // devstd: sqrt( (⟨x²⟩ − ⟨x⟩²) / (i−1) )
        devstd.push_back(sqrt((SquaredCumSum / i - pow(mean.at(i-1), 2)) / (i-1)));

        // debug: confronto SquaredCumSum/i e pow(mean.at(i-1),2)
        // cout << endl << SquaredCumSum/i << "   " << pow(mean.at(i-1),2) << endl;
    }

    // Salvo su file "data.txt"
    outData(mean, devstd, "data.txt");
}

// Calcola la varianza di un singolo blocco di 'throws' estrazioni
double OneBlockVariance(Random* rnd, unsigned int throws){
    double x;                     // variabile per l'estrazione
    double cum_sum = 0;           // somma cumulativa dei valori
    double squared_cum_sum = 0;   // somma cumulativa dei quadrati
    double variance;              // varianza del blocco

    for (unsigned int i = 0; i < throws; i++) {
        x = rnd->Rannyu();        // estraggo un numero casuale
        cum_sum += x;
        squared_cum_sum += pow(x, 2);
    }
    // varianza = ⟨x²⟩ − ⟨x⟩²
    variance = (squared_cum_sum / throws) - pow((cum_sum / throws), 2);
    return variance;             // restituisco la varianza
}

// Costruisce progressivamente la media della varianza e la sua dev. std. per blocchi
void VarianceAndDevStd(Random* rnd, unsigned int nBlocksMax, unsigned int nThrows){
    double CumSum = 0;             
    double SquaredCumSum = 0;      
    double x;                      
    vector<double> mean, devstd;   // media della varianza e dev. std.

    // Primo blocco: devstd = 0
    x = OneBlockVariance(rnd, nThrows);
    CumSum += x;
    SquaredCumSum += pow(x, 2);
    mean.push_back(CumSum);
    devstd.push_back(0);

    // Blocchi successivi
    for (unsigned int i = 2; i <= nBlocksMax; i++) {
        x = OneBlockVariance(rnd, nThrows);
        CumSum += x;
        SquaredCumSum += pow(x, 2);
        mean.push_back(CumSum / double(i));
        devstd.push_back(sqrt((SquaredCumSum / i - pow(mean.at(i-1), 2)) / (i-1)));

        // debug:
        // cout << endl << SquaredCumSum/i << "   " << pow(mean.at(i-1),2) << endl;
    }

    // Salvo su file "dataSigma.txt"
    outData(mean, devstd, "dataSigma.txt");
}

// Conta quante estrazioni finiscono in ciascuno dei 100 bin uniformi [0,1)
vector<unsigned int> Vector_Uniform_Distrib(Random* rnd, unsigned int N){
    double x;                        
    unsigned int index;             // indice del bin (0–99)
    vector<unsigned int> counter(100, 0); // inizializzo 100 contatori a zero

    for (unsigned int i = 0; i < N; i++) {
        x = rnd->Rannyu();           // numero uniforme [0,1)
        index = int(x * 100);        // bin = floor(x*M)
        counter.at(index) += 1;      // incremento il conteggio
    }
    return counter;                 
}

// Calcola il chi quadro: Σ (osservato − atteso)² / atteso
double chi_squared(vector<unsigned int> counter, unsigned int nThrows){
    double chi2 = 0;
    double nM = nThrows / counter.size(); // atteso per bin = nThrows/M

    for (unsigned int i = 0; i < counter.size(); i++) {
        chi2 += pow((counter.at(i) - nM), 2) / nM;
    }
    return chi2;
}

// Ripete nTot volte il test del chi quadro e salva i risultati divisi per nTot
void outChi(Random* rnd, unsigned int nTot, unsigned int nThrows, const char file[]){
    vector<unsigned int> vCounter; 
    vector<double> ChiSquared;      

    // Simulazione chi quadro nTot volte
    for (unsigned int i = 0; i < nTot; i++) {
        vCounter = Vector_Uniform_Distrib(rnd, nThrows);
        ChiSquared.push_back(chi_squared(vCounter, nThrows));
    }

    // Scrivo su file
    ofstream output;
    output.open(file);
    if (output.fail()) {
        cout << endl
             << "Errore di apertura dello stream di output, "
                "il file su cui vuoi trascrivere è nella stessa cartella dell'eseguibile?"
             << endl;
        return;
    }
    for (unsigned int i = 0; i < ChiSquared.size(); i++) {
        // divido il valore per nTot per normalizzare
        output << ChiSquared.at(i) / nTot << endl;
    }
    output.close();
}
