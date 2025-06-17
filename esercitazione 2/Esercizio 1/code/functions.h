// functions.h

#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <vector>

// Campionamento uniforme su [a,b]: calcola ⟨f⟩ per blocco di nX punti
double OneBlockMean(Random* rnd, unsigned int nX, double a, double b);

// Importanza con p(x)=2(1–x): calcola ⟨f/p⟩ per blocco di nX punti
double OneBlockMeanPx(Random* rnd, unsigned int nX);

// Stampa mean e StdDev su file in due colonne
void outData(std::vector<double> mean, std::vector<double> StdDev, const char file[]);

// Calcola media e varianza dei blocchi con campionamento uniforme
void MeanAndVariance(Random* rnd, unsigned int nBlocks, unsigned int nX, double a, double b);

// Calcola media e varianza dei blocchi con importanza (p(x)≈1−x)
void MeanAndVariancePx(Random* rnd, unsigned int nBlocks, unsigned int nX);

#endif
