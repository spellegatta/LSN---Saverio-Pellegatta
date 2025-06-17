// main.cpp

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "functions.h"

using namespace std;

int main (int argc, char *argv[]){

   Random rnd;                 // istanza del generatore
   int seed[4];
   int p1, p2;

   // Lettura due numeri primi da "Primes"
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   // Lettura del seme da "seed.in"
   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while (!input.eof()) {
         input >> property;
         if (property == "RANDOMSEED") {
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed, p1, p2);  // inizializza il generatore
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   // Simulazione: 100 blocchi da 1e4 punti ciascuno
   MeanAndVariance(&rnd, 100, 1E4, 0, 1);
   MeanAndVariancePx(&rnd, 100, 1E4);

   rnd.SaveSeed();            // salvo lo stato del generatore
   return 0;
}
