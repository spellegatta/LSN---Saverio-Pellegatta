// functions.cpp

#include <vector>
#include <iostream>
#include "random.h"
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

// Scrive su file i vettori mean (colonna 1) e StdDev (colonna 2)
void outData(std::vector<double> mean, std::vector<double> StdDev, const char file[]){
    ofstream output;
    output.open(file);
    if(output.fail()){
        cout << endl
             << "Errore di apertura dello stream di output, "
                "il file su cui vuoi trascrivere è nella stessa cartella dell'eseguibile?"
             << endl;
        return;                   // esco in caso di errore
    }
    for(unsigned int i = 0; i < mean.size(); i++){
        output << mean.at(i)
               << setw(12)       // spaziatura per allineare StdDev
               << StdDev.at(i)
               << endl;
    }
    output.close();             // chiudo lo stream
}

// Esegue un singolo esperimento di Buffon e restituisce stima di π
double experiment(Random* rnd, unsigned int nThrows, double L, double d) {
    if (L > d) {
        // Riformulato: controllo preliminare sui parametri
        std::cerr << "Errore: L (" << L 
                  << ") deve essere <= distanza tra le linee d (" 
                  << d << ")\n";
        return 0.0;
    }

    unsigned int HitCounter = 0;  // conta le intersezioni ago-linee

    for (unsigned int i = 0; i < nThrows; ++i) {
        // 1) distanza del centro ago dalla linea più vicina [0, d/2)
        double x = rnd->Rannyu(0.0, d/2.0);

        // 2) generazione di un angolo θ tramite punto (a,b) nel cerchio unitario
        double a, b, r;
        do {
            a = rnd->Rannyu(-1.0, 1.0);
            b = rnd->Rannyu(-1.0, 1.0);
            r = a*a + b*b;
        } while (r > 1.0 or r == 0.0);
        // r = a² + b² ∈ (0,1] assicura una direzione uniforme

        // 3) sin(theta) = componente y normalizzata su r
        double sin_theta = fabs(b) / sqrt(r);

        // 4) condizione di intersezione: x <= (L/2)·sin(theta)
        if (x <= (L / 2.0) * sin_theta) {
            HitCounter++;
        }
    }

    if (HitCounter == 0) {
        std::cerr << "Nessun colpo: impossibile stimare π\n";
        return 0.0;
    }

    // 5) formula di Buffon: π ≈ (2·L·N) / (Hits·d)
    return (2.0 * L * nThrows) / (HitCounter * d);
}

// Calcola media e deviazione standard della stima di π su blocchi
void MeanAndDevStd(Random* rnd, unsigned int nBlocksMax, unsigned int nThrows, double L, double d){
    double CumSum = 0;             // somma cumulativa delle stime
    double SquaredCumSum = 0;      // somma cumulativa dei quadrati
    double x;                      // variabile di appoggio
    vector<double> mean, devstd;   // memorizzano media e dev.std per i blocchi

    // Primo blocco: devstd = 0
    x = experiment(rnd, nThrows, L, d);
    CumSum += x;
    SquaredCumSum += pow(x, 2);
    mean.push_back(CumSum);        // media cumulativa con 1 blocco
    devstd.push_back(0);           // deviazione standard nulla

    // Blocchi successivi da 2 a nBlocksMax
    for (unsigned int i = 2; i <= nBlocksMax; i++){
        x = experiment(rnd, nThrows, L, d);
        CumSum += x;
        SquaredCumSum += pow(x, 2);
        mean.push_back(CumSum / double(i));
        // devstd: sqrt((⟨x²⟩ − ⟨x⟩²)/(i−1))
        devstd.push_back(sqrt((SquaredCumSum / i - pow(mean.at(i-1), 2)) / (i-1)));

        // debug:
        // cout << SquaredCumSum/i << "   " << pow(mean.at(i-1),2) << endl;
    }

    // Salvo risultati su file "data.txt"
    outData(mean, devstd, "data.txt");
}
