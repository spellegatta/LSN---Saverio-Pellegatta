// main.cpp

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "functions.h"

using namespace std;

int main (int argc, char *argv[]){

   Random rnd;                     // istanza del generatore Mersenne-Twister
   int seed[4];
   int p1, p2;

   // Lettura dei primi per l'inizializzazione
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
            rnd.SetRandom(seed, p1, p2);  // setto il generatore
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   // Eseguo la simulazione di Buffon: 1000 blocchi da 1e6 estrazioni
   MeanAndDevStd(&rnd, 1000, 1E6, 0.001, 0.1);

   rnd.SaveSeed();                // salvo lo stato del generatore
   return 0;
}
