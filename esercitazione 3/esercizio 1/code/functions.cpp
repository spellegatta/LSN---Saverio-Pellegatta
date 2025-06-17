// functions.cpp

#include <vector>
#include <iostream>
#include "random.h"
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

// Simula S(t) in un solo colpo: S₀·exp((μ−σ²/2)t + σ·Wₜ), con Wₜ ~ N(0,t)
double St_function(Random* rnd, double S0, double mu, double sigma, double t){
    double W = rnd->Gauss(0, sqrt(t));  
    return S0 * exp((mu - pow(sigma, 2) / 2.) * t + sigma * W);
}

// Passo incrementale di S_t→S_{t+Δt}: S·exp((μ−σ²/2)Δt + σ·Z√Δt), Z~N(0,1)
double St_function_step(Random* rnd, double S_t_i, double mu, double sigma, double DeltaT){
    double Z = rnd->Gauss(0, 1);
    return S_t_i * exp((mu - pow(sigma, 2) / 2.) * DeltaT + sigma * Z * sqrt(DeltaT));
}

// Itera dal tempo t₀ a T in nStep passi, partendo da S₀, usando St_function_step
double Iterate_S_t_i(Random* rnd, double nStep, double t0, double T, double mu, double sigma, double S0){
    double S_i1, S_i2;
    double DeltaT = (T - t0) / nStep;          // ampiezza di ciascun passo temporale
    S_i1 = St_function(rnd, S0, mu, sigma, DeltaT);  // primo passo da S₀

    for (unsigned int i = 0; i < (nStep - 1); i++){
        S_i2 = St_function_step(rnd, S_i1, mu, sigma, DeltaT);
        S_i1 = S_i2;                           // aggiorno per il passo successivo
    }
    return S_i1;                              // restituisco S_T simulato
}

// Payoff di una call europea: e^{−μT}·max(0, S_T − K)
double Call(double mu, double T, double K, double ST){
    return exp(-mu * T) * max(0., ST - K);
}

// Payoff di una put europea: e^{−μT}·max(0, K − S_T)
double Put(double mu, double T, double K, double ST){
    return exp(-mu * T) * max(0., K - ST);
}

// Scrive su file mean e StdDev in due colonne
void outData(std::vector<double> mean, std::vector<double> StdDev, const char file[]){
    ofstream output;
    output.open(file);
    if (output.fail()){
        cout << endl
             << "Errore di apertura dello stream di output, "
                "il file su cui vuoi trascrivere è nella stessa cartella dell'eseguibile?"
             << endl;
        return;                              // esco se non apro il file
    }
    for (unsigned int i = 0; i < mean.size(); i++){
        output << mean.at(i)
               << setw(12)                   // spaziatura colonna StdDev
               << StdDev.at(i)
               << endl;
    }
    output.close();                          // chiudo lo stream
}

// Calcola su un blocco la media di call e put, ritorna vettore [⟨call⟩, ⟨put⟩]
vector<double> OneBlockMean(Random* rnd, unsigned int throws,
                            double S0, double mu, double sigma, double t, double K){
    double x, y, ST;                          // variabili d’appoggio
    double CumSumX = 0, CumSumY = 0;          // somme cumulative call e put
    vector<double> results(2);

    for (unsigned int i = 0; i < throws; i++){
        ST = St_function(rnd, S0, mu, sigma, t);  // simulo S_T
        x = Call(mu, t, K, ST);                   // payoff call
        y = Put(mu, t, K, ST);                    // payoff put
        CumSumX += x;                             // accumulo call
        CumSumY += y;                             // accumulo put
    }
    results.at(0) = CumSumX / throws;            // media call
    results.at(1) = CumSumY / throws;            // media put

    return results;
}

// Costruisce progressivamente media e dev.std. di call e put su nBlocksMax blocchi
void CallAndPut(Random* rnd, unsigned int nBlocksMax, unsigned int nThrows,
                double S0, double mu, double sigma, double t, double K){
    double CumSumCall = 0, SquaredCumSumCall = 0;
    double CumSumPut  = 0, SquaredCumSumPut  = 0;
    vector<double> call, devstd_call, put, devstd_put;
    vector<double> x;                         // vettore [call, put] di un blocco

    // Primo blocco: devstd = 0
    x = OneBlockMean(rnd, nThrows, S0, mu, sigma, t, K);
    CumSumCall += x.at(0);
    SquaredCumSumCall += pow(x.at(0), 2);
    CumSumPut += x.at(1);
    SquaredCumSumPut += pow(x.at(1), 2);
    call.push_back(CumSumCall);
    devstd_call.push_back(0);
    put.push_back(CumSumPut);
    devstd_put.push_back(0);

    // Blocchi successivi da 2 a nBlocksMax
    for (unsigned int i = 2; i <= nBlocksMax; i++){
        x = OneBlockMean(rnd, nThrows, S0, mu, sigma, t, K);
        CumSumCall += x.at(0);
        SquaredCumSumCall += pow(x.at(0), 2);
        CumSumPut += x.at(1);
        SquaredCumSumPut += pow(x.at(1), 2);

        call.push_back(CumSumCall / double(i));
        devstd_call.push_back(sqrt((SquaredCumSumCall / i
                                    - pow(call.at(i-1), 2))
                                   / (i-1)));

        put.push_back(CumSumPut / double(i));
        devstd_put.push_back(sqrt((SquaredCumSumPut / i
                                   - pow(put.at(i-1), 2))
                                  / (i-1)));
    }

    // Salvo i risultati
    outData(call,     devstd_call, "dataCall.txt");
    outData(put,      devstd_put,  "dataPut.txt");
}

// Simile a OneBlockMean ma usando Iterate_S_t_i per cammini temporali a step
vector<double> OneBlockMean_iter(Random* rnd, unsigned int throws,
                                 double S0, double mu, double sigma,
                                 double t0, double T, double K, double nStep){
    double x, y, ST;
    double CumSumX = 0, CumSumY = 0;
    vector<double> results(2);

    for (unsigned int i = 0; i < throws; i++){
        ST = Iterate_S_t_i(rnd, nStep, t0, T, mu, sigma, S0);
        x = Call(mu, T, K, ST);
        y = Put(mu, T, K, ST);
        CumSumX += x;
        CumSumY += y;
    }
    results.at(0) = CumSumX / throws;
    results.at(1) = CumSumY / throws;
    return results;
}

// Costruisce media e dev.std. di call/put usando Iterate_S_t_i
void CallAndPut_iter(Random* rnd, unsigned int nBlocksMax, unsigned int nThrows,
                     double S0, double mu, double sigma, double t0,
                     double T, double K, unsigned int nStep){
    double CumSumCall = 0, SquaredCumSumCall = 0;
    double CumSumPut  = 0, SquaredCumSumPut  = 0;
    vector<double> call, devstd_call, put, devstd_put;
    vector<double> x;

    // Primo blocco
    x = OneBlockMean_iter(rnd, nThrows, S0, mu, sigma, t0, T, K, nStep);
    CumSumCall += x.at(0);
    SquaredCumSumCall += pow(x.at(0), 2);
    CumSumPut += x.at(1);
    SquaredCumSumPut += pow(x.at(1), 2);
    call.push_back(CumSumCall);
    devstd_call.push_back(0);
    put.push_back(CumSumPut);
    devstd_put.push_back(0);

    // Blocchi successivi
    for (unsigned int i = 2; i <= nBlocksMax; i++){
        x = OneBlockMean_iter(rnd, nThrows, S0, mu, sigma, t0, T, K, nStep);
        CumSumCall += x.at(0);
        SquaredCumSumCall += pow(x.at(0), 2);
        CumSumPut += x.at(1);
        SquaredCumSumPut += pow(x.at(1), 2);

        call.push_back(CumSumCall / double(i));
        devstd_call.push_back(sqrt((SquaredCumSumCall / i
                                    - pow(call.at(i-1), 2))
                                   / (i-1)));

        put.push_back(CumSumPut / double(i));
        devstd_put.push_back(sqrt((SquaredCumSumPut / i
                                   - pow(put.at(i-1), 2))
                                  / (i-1)));
    }

    // Salvo i risultati iterativi
    outData(call,     devstd_call,  "dataCall_iter.txt");
    outData(put,      devstd_put,   "dataPut_iter.txt");
}
