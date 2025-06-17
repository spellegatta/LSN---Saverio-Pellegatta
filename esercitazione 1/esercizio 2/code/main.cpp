// main.cpp

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "functions.h"

using namespace std;

int main (int argc, char *argv[]){

   Random rnd;                       // istanza del generatore
   int seed[4];
   int p1, p2;

   // Legge due numeri primi da file "Primes"
   ifstream Primes("Primes");
   if (Primes.is_open()) {
      Primes >> p1 >> p2;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   // Legge il seme da file "seed.in"
   ifstream input("seed.in");
   string property;
   if (input.is_open()) {
      while (!input.eof()) {
         input >> property;
         if (property == "RANDOMSEED") {
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed, p1, p2);  // inizializza il generatore
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   // Vettore con i diversi numeri di estrazioni da testare
   vector<unsigned int> nThrows = {1, 2, 10, 100};

   // Per ogni valore di nThrows, genero 1e4 medie e le salvo su file
   for (unsigned int i = 0; i < nThrows.size(); i++){
      vector<double> vUniformMean, vExponentialMean, vCLMean;

      for (unsigned int j = 0; j < 1E4; j++) {
         vUniformMean.push_back(   UniformMean(&rnd, nThrows.at(i)) );
         vExponentialMean.push_back( ExponentialMean(&rnd, nThrows.at(i), 1) );
         vCLMean.push_back(        CauchyLorentzMean(&rnd, nThrows.at(i), 0, 1) );
      }

      // Costruisco i nomi dei file in output e salvo i dati
      string filenameUnif = "dataUniform"   + to_string(nThrows.at(i)) + ".txt";
      outData(vUniformMean, filenameUnif.c_str());

      string filenameExp  = "dataExponential" + to_string(nThrows.at(i)) + ".txt";
      outData(vExponentialMean, filenameExp.c_str());

      string filenameCL   = "dataCL"          + to_string(nThrows.at(i)) + ".txt";
      outData(vCLMean, filenameCL.c_str());
   }

   rnd.SaveSeed();                   // salva lo stato del generatore
   return 0;
}
