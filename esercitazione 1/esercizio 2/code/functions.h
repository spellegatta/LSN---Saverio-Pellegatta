// functions.h

#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <vector>

// Calcola e restituisce la media di numeri generati casualmente con distribuzione uniforme
double UniformMean(Random* rnd, unsigned int nThrows);

// Calcola e restituisce la media di numeri generati con distribuzione esponenziale (parametro lambda)
double ExponentialMean(Random* rnd, unsigned int nThrows, double lambda);

// Calcola e restituisce la media di numeri generati con distribuzione Cauchy-Lorentz (parametri mu, gamma)
double CauchyLorentzMean(Random* rnd, unsigned int nThrows, double mu, double gamma);

// Scrive su file i valori contenuti in vec, uno per riga
void outData(std::vector<double> vec, const char file[]);

#endif
