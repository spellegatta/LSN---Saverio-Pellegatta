// functions.cpp

#include <vector>
#include <iostream>
#include "random.h"
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

// Calcola la media di un blocco di nX punti per l'integrazione di f(x) su [a,b]
// mediante campionamento uniforme x = a + (b–a)*Rannyu(), e restituisce ⟨f⟩
double OneBlockMean(Random* rnd, unsigned int nX, double a, double b){
    double fx, x;                  // fx = f(x), x = punto campionato
    double cum_sum = 0;            // somma cumulativa dei valori f(x)
    for (unsigned int i = 0; i < nX; i++){
        x = rnd->Rannyu() * (b - a) + a;            // campiono x ∈ [a,b)
        fx = (M_PI/2.) * cos(M_PI * x / 2.);        // f(x) = (π/2)·cos(πx/2)
        cum_sum += fx;                              // accumulo la somma
    }
    return cum_sum / nX;            // ritorno la media del blocco
}

// Calcola la media di un blocco di nX punti usando importanza con p(x)=2(1–x)
// campionando x = 1 – sqrt(1 – Rannyu()), e restituisce ⟨f(x)/p(x)⟩
double OneBlockMeanPx(Random* rnd, unsigned int nX){
    double fx, x;                  // fx = peso*f(x), x = punto campionato
    double cum_sum = 0;            // somma cumulativa dei valori fx
    for (unsigned int i = 0; i < nX; i++){
        x = 1 - sqrt(1 - rnd->Rannyu());            // metodo inverso per p(x)
        fx = (M_PI/4.) * cos(M_PI * x / 2.) / (1 - x);
        cum_sum += fx;                              // accumulo la somma
    }
    return cum_sum / nX;            // ritorno la media del blocco
}

// Scrive su file i vettori mean (colonna 1) e StdDev (colonna 2)
void outData(std::vector<double> mean, std::vector<double> StdDev, const char file[]){
    ofstream output;
    output.open(file);
    if (output.fail()){
        cout << endl
             << "Errore di apertura dello stream di output, "
                "il file su cui vuoi trascrivere è nella stessa cartella dell'eseguibile?"
             << endl;
        return;                   // esco in caso di errore di apertura
    }
    for (unsigned int i = 0; i < mean.size(); i++){
        output << mean.at(i)
               << setw(12)       // allineo la seconda colonna
               << StdDev.at(i)
               << endl;
    }
    output.close();             // chiudo lo stream
}

// Calcola media e deviazione standard dei blocchi per l'integrazione uniforme
void MeanAndVariance(Random* rnd, unsigned int nBlocks, unsigned int nX, double a, double b){
    double CumSum = 0;           // somma cumulativa delle medie
    double SquaredCumSum = 0;    // somma cumulativa dei quadrati delle medie
    double x;                    // variabile di supporto
    vector<double> mean, devstd; // vettori di output: media cumulativa e dev.std.

    // Assicuro a <= b; se no scambio
    if (b < a){
        double c = a;
        a = b;
        b = c;
    }

    // Primo blocco: devstd = 0
    x = (b - a) * OneBlockMean(rnd, nX, a, b);
    CumSum += x;
    SquaredCumSum += pow(x, 2);
    mean.push_back(CumSum);      // media cumulativa con 1 blocco
    devstd.push_back(0);         // deviazione standard nulla

    // Blocchi successivi da 2 a nBlocks
    for (unsigned int i = 2; i <= nBlocks; i++){
        x = (b - a) * OneBlockMean(rnd, nX, a, b);
        CumSum += x;
        SquaredCumSum += pow(x, 2);
        mean.push_back(CumSum / double(i));
        // devstd = (b–a)·sqrt((⟨x²⟩−⟨x⟩²)/(i−1))
        devstd.push_back((b - a) * sqrt((SquaredCumSum / i - pow(mean.at(i-1), 2)) / (i-1)));
    }

    // Salvo su file "data.txt"
    outData(mean, devstd, "data.txt");
}

// Calcola media e deviazione standard dei blocchi per l'integrazione importance-sampling
void MeanAndVariancePx(Random* rnd, unsigned int nBlocks, unsigned int nX){
    double CumSum = 0;           // somma cumulativa delle medie
    double SquaredCumSum = 0;    // somma cumulativa dei quadrati delle medie
    double x;                    // variabile di supporto
    vector<double> mean, devstd; // vettori di output

    // Primo blocco: devstd = 0
    x = OneBlockMeanPx(rnd, nX);
    CumSum += x;
    SquaredCumSum += pow(x, 2);
    mean.push_back(CumSum);
    devstd.push_back(0);

    // Blocchi successivi
    for (unsigned int i = 2; i <= nBlocks; i++){
        x = OneBlockMeanPx(rnd, nX);
        CumSum += x;
        SquaredCumSum += pow(x, 2);
        mean.push_back(CumSum / double(i));
        devstd.push_back(sqrt((SquaredCumSum / i - pow(mean.at(i-1), 2)) / (i-1)));
    }

    // Salvo su file "dataPx.txt"
    outData(mean, devstd, "dataPx.txt");
}
