// main.cpp

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "functions.h"
#include "RandomWalk.h"

using namespace std;

int main (int argc, char *argv[]){

   Random rnd;                          // generatore di numeri casuali
   int seed[4];
   int p1, p2;

   // Lettura dei primi da file "Primes"
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   // Lettura del seme da "seed.in"
   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while (!input.eof()){
         input >> property;
         if (property == "RANDOMSEED"){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed, p1, p2); // inizializzo il generatore
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   // Simulazione RW reticolare
   RandomWalk* rw = new RandomWalk(&rnd);
   MeanAndDevStd(rw, 100, 1000, 100, "ReticularRandomWalk.txt");
   delete rw;

   // Simulazione RW continuo
   RandomWalk* rw_continuum = new RandomWalk_Continuum(&rnd);
   MeanAndDevStd(rw_continuum, 100, 1000, 100, "ContinousRandomWalk.txt");
   delete rw_continuum;

   rnd.SaveSeed();                     // salvo stato del generatore
   return 0;
}
