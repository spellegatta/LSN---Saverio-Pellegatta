#include <vector>
#include <iostream>
#include "random.h"
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

double OneBlockMean(Random* rnd, unsigned int throws){
    double x; //variabile di appoggio
    double cum_sum=0; //somma cumulativa
    for (unsigned int i=0; i<throws; i++){
        x=rnd->Rannyu(); //estraggo un numero casuale
        cum_sum+=x; //somma cumulativa
    }
    return cum_sum/throws; //ritorna la media sul singolo blocco
}

void outData(std::vector<double> mean, std::vector<double> StdDev, const char file[]){
    ofstream output;
    output.open(file);
    if(output.fail()){
        cout << endl << "Errore di apertura dello stream di output, il file su cui vuoi trascrivere è nella stessa cartella dell'eseguibile?" << endl;
        return;
    }
    for(unsigned int i=0; i<mean.size(); i++){
        output << mean.at(i) << setw(12) << StdDev.at(i) << endl;
    }
    output.close();

}

void MeanAndDevStd(Random* rnd, unsigned int nBlocksMax, unsigned int nThrows){
    double CumSum=0; //somma cumulativa delle medie dei blocchi
    double SquaredCumSum=0; //somma cumulativa dei quadrati delle medie dei blocchi
    double x; //variabile di appoggio
    vector<double> mean, devstd; //vettori dove scrivo le medie e le dev standard dei gruppi di blocchi

    //se ho solo un blocco impongo la devstd nulla
    x=OneBlockMean(rnd, nThrows);
    CumSum+=x;
    SquaredCumSum+=pow(x,2);
    mean.push_back(CumSum);
    devstd.push_back(0);
    
    //calcolo media e dev std con un numero di blocchi maggiore di 1
    for (unsigned int i=2; i<=nBlocksMax; i++){
        x=OneBlockMean(rnd, nThrows);
        CumSum+=x;
        SquaredCumSum+=pow(x,2);
        mean.push_back(CumSum/double(i));
        devstd.push_back(sqrt((SquaredCumSum/i-pow(mean.at(i-1),2))/(i-1)));

        //cout << endl << SquaredCumSum/i << "   " << pow(mean.at(i-1),2) << endl;


    }

    //Stampo i dati su file
    outData(mean,devstd,"data.txt");

    
}


double OneBlockVariance(Random* rnd, unsigned int throws){
    double x; //variabile di appoggio
    double cum_sum=0; //somma cumulativa
    double squared_cum_sum=0; //somma cumulativa dei quadrati
    double variance; //varianza
    for (unsigned int i=0; i<throws; i++){
        x=rnd->Rannyu(); //estraggo un numero casuale
        cum_sum+=x; //somma cumulativa
        squared_cum_sum+=pow(x,2); // somma cumulativa dei quadrati
    }
    variance=(squared_cum_sum/throws)-pow((cum_sum/throws),2);
    return variance; //ritorna la varianza
}

void VarianceAndDevStd(Random* rnd, unsigned int nBlocksMax, unsigned int nThrows){
    double CumSum=0; //somma cumulativa delle medie dei blocchi
    double SquaredCumSum=0; //somma cumulativa dei quadrati delle medie dei blocchi
    double x; //variabile di appoggio
    vector<double> mean, devstd; //vettori dove scrivo le medie e le dev standard dei gruppi di blocchi

    //se ho solo un blocco impongo la devstd nulla
    x=OneBlockVariance(rnd, nThrows);
    CumSum+=x;
    SquaredCumSum+=pow(x,2);
    mean.push_back(CumSum);
    devstd.push_back(0);
    
    //calcolo media e dev std con un numero di blocchi maggiore di 1
    for (unsigned int i=2; i<=nBlocksMax; i++){
        x=OneBlockVariance(rnd, nThrows);
        CumSum+=x;
        SquaredCumSum+=pow(x,2);
        mean.push_back(CumSum/double(i));
        devstd.push_back(sqrt((SquaredCumSum/i-pow(mean.at(i-1),2))/(i-1)));

        //cout << endl << SquaredCumSum/i << "   " << pow(mean.at(i-1),2) << endl;


    }

    //Stampo i dati su file
    outData(mean,devstd,"dataSigma.txt");

    
}

vector<unsigned int> Vector_Uniform_Distrib(Random* rnd, unsigned int N){
    double x;//variabile di appoggio
    unsigned int index; //indice 
    vector<unsigned int> counter(100, 0); //inizializzo il counter con 100 posti a 0
    for(unsigned int i = 0; i<N; i++){
        x=rnd->Rannyu(); //estraggo un numero tra 0 e 1 con distribuzione uniforme
        index=int(x*100); // se moltiplico per il numero di bins il numero e lo arrotondo verso il basso ottengo l'indice del vettore da incrementare
        counter.at(index)+=1; //incremento il counter
    }
    return counter;
}
double chi_squared(vector<unsigned int> counter, unsigned int nThrows){
    double chi_squared=0;
    double nM=nThrows/counter.size();// n/M
    for(unsigned int i=0; i<counter.size();i++){
        chi_squared+=pow((counter.at(i)-nM),2)/nM;
    }
    return chi_squared;
}

void outChi(Random* rnd, unsigned int nTot, unsigned int nThrows, const char file[]){
    vector<unsigned int> vCounter; //vettore che conta le estrazioni
    vector<double> ChiSquared; //chi quadro

    //simulo
    for (unsigned int i=0; i<nTot; i++){
        vCounter=Vector_Uniform_Distrib(rnd, nThrows);
        ChiSquared.push_back(chi_squared(vCounter, nThrows));
    }

    //stampo su file
    ofstream output;
    output.open(file);
    if(output.fail()){
        cout << endl << "Errore di apertura dello stream di output, il file su cui vuoi trascrivere è nella stessa cartella dell'eseguibile?" << endl;
        return;
    }
    for(unsigned int i=0; i<ChiSquared.size(); i++){
        output << ChiSquared.at(i)<< endl;
    }
    output.close();

}




