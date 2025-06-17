// functions.h

#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <vector>

// Esegue un esperimento di Buffon e restituisce la stima di π
double experiment(Random* rnd, unsigned int nThrows, double L, double d);

// Calcola media e dev.std. della stima di π, su nBlocksMax blocchi
void MeanAndDevStd(Random* rnd, unsigned int nBlocksMax,
                   unsigned int nThrows, double L, double d);

// Scrive mean e StdDev su file in due colonne
void outData(std::vector<double> mean, std::vector<double> StdDev,
             const char file[]);

#endif
