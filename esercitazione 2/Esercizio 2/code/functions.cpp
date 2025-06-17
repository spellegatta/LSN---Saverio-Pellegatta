// functions.cpp

#include <vector>
#include <iostream>
#include "random.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include "RandomWalk.h"

using namespace std;

// Calcola la media dei vettori di distanza dal punto di partenza
// su un singolo blocco: lancia 'throws' cammini di 'nSteps' passi.
vector<double> OneBlockMean(RandomWalk* rw, unsigned int throws, unsigned int nSteps){
    vector<double> x;                        // vettore di posizioni per il singolo cammino
    vector<double> cum_sum(nSteps, 0);       // somma cumulativa delle posizioni per ogni passo

    for (unsigned int i = 0; i < throws; i++){
        x = rw->GenerateWalk(nSteps);        // genero un cammino di nSteps passi
        for (unsigned int j = 0; j < x.size(); j++){
            cum_sum.at(j) += x.at(j);        // accumulo la distanza al passo j
        }
    }
    // Divido per il numero di lanci per ottenere la media di ciascun passo
    for (unsigned int j = 0; j < cum_sum.size(); j++){
        cum_sum.at(j) = cum_sum.at(j) / throws;
    }
    return cum_sum;                          // restituisco il vettore delle medie
}

// Scrive su file mean (colonna 1) e StdDev (colonna 2) per ogni passo
void outData(std::vector<double> mean, std::vector<double> StdDev, const char file[]){
    ofstream output;
    output.open(file);
    if (output.fail()){
        cout << endl
             << "Errore di apertura dello stream di output, "
                "il file su cui vuoi trascrivere è nella stessa cartella dell'eseguibile?"
             << endl;
        return;                            // esco se non riesco ad aprire il file
    }
    for (unsigned int i = 0; i < mean.size(); i++){
        output << mean.at(i)
               << setw(12)               // spaziatura per la colonna StdDev
               << StdDev.at(i)
               << endl;
    }
    output.close();                        // chiudo lo stream
}

// Calcola media e deviazione standard progressiva sui blocchi
// per ogni passo del cammino e salva i risultati su 'filename'
void MeanAndDevStd(RandomWalk* rw, unsigned int nBlocksMax, unsigned int nThrows, unsigned int nSteps, const char* filename){
    vector<double> CumSum(nSteps, 0);       // somma cumulativa delle medie per passo
    vector<double> SquaredCumSum(nSteps, 0);// somma cumulativa dei quadrati delle medie
    vector<double> x;                       // vettore temporaneo per il blocco
    vector<double> mean(nSteps, 0);         // vettore delle medie finali
    vector<double> devstd(nSteps, 0);       // vettore delle dev. std. finali

    // Accumulo su tutti i blocchi
    for (unsigned int i = 1; i <= nBlocksMax; i++){
        x = OneBlockMean(rw, nThrows, nSteps);
        for (unsigned int j = 0; j < x.size(); j++){
            CumSum.at(j)       += x.at(j);       // accumulo media
            SquaredCumSum.at(j) += pow(x.at(j), 2); // accumulo quadrato media
        }
    }

    // Calcolo media e dev.std. per ciascun passo
    for (unsigned int j = 0; j < SquaredCumSum.size(); j++){
        mean.at(j)   = CumSum.at(j) / double(nBlocksMax);
        devstd.at(j) = sqrt(
            (SquaredCumSum.at(j) / double(nBlocksMax)
             - pow(CumSum.at(j) / double(nBlocksMax), 2))
            / (nBlocksMax - 1)
        );
    }

    // Salvo su file
    outData(mean, devstd, filename);
}
