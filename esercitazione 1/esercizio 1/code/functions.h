// functions.h

#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <vector>

double OneBlockMean(Random* rnd, unsigned int throws);  // media di un blocco: A_i
void MeanAndDevStd(Random* rnd, unsigned int nBlocksMax, unsigned int nThrows);  
    // media e dev.std progressive, output su data.txt
void outData(std::vector<double> mean, std::vector<double> StdDev, const char file[]);
    // stampa mean e stddev su file in due colonne
double OneBlockVariance(Random* rnd, unsigned int throws);  
    // varianza di un blocco
void VarianceAndDevStd(Random* rnd, unsigned int nBlocksMax, unsigned int nThrows);
    // media della varianza e dev.std, output su dataSigma.txt
std::vector<unsigned int> Vector_Uniform_Distrib(Random* rnd, unsigned int N);
    // conta N estrazioni uniformi in 100 bin
double chi_squared(std::vector<unsigned int> counter, unsigned int nThrows);
    // calcolo chi quadro con counter e numero estrazioni
void outChi(Random* rnd, unsigned int nTot, unsigned int nThrows, const char file[]);
    // simula outChi nTot volte, salva in file

#endif
