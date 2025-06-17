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

   // Leggo i due numeri primi da file "Primes"
   ifstream Primes("Primes");
   if (Primes.is_open()) {
      Primes >> p1 >> p2;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   // Leggo il seme da "seed.in"
   ifstream input("seed.in");
   string property;
   if (input.is_open()) {
      while (!input.eof()) {
         input >> property;
         if (property == "RANDOMSEED") {
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed, p1, p2);  // inizializzo il generatore
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   // Simulazioni principali:
   MeanAndDevStd(&rnd, 500, 1e5);      // media e devstd su 500 blocchi di 1e5 estrazioni
   VarianceAndDevStd(&rnd, 500, 1e5);  // varianza e devstd su 500 blocchi
   outChi(&rnd, 100, 1e5, "dataChi.txt"); // chi² su 100 ripetizioni

   rnd.SaveSeed();              // salvo lo stato del generatore
   return 0;
}

