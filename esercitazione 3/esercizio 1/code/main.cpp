// main.cpp

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "functions.h"

using namespace std;

int main (int argc, char *argv[]){

   Random rnd;                // istanza del generatore
   int seed[4];
   int p1, p2;

   // Lettura dei primi
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   // Lettura del seme
   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while (!input.eof()){
         input >> property;
         if (property == "RANDOMSEED"){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed, p1, p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   // Parametri dell’opzione
   double T     = 1;
   double S0    = 100;
   double K     = 100;
   double mu    = 0.1;
   double sigma = 0.25;
   double t0    = 0;
   unsigned int nStep = 100;

   // Simulazione europea (senza discretizzazione)
   CallAndPut(&rnd,       100, 1E6, S0, mu, sigma, T, K);

   // Simulazione con discretizzazione temporale
   CallAndPut_iter(&rnd,  100, 1E6, S0, mu, sigma, t0, T, K, nStep);

   rnd.SaveSeed();         // salvo lo stato del generatore
   return 0;
}
