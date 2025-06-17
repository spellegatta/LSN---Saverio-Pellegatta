// functions.h

#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <vector>
#include "RandomWalk.h"

// Genera 'throws' cammini di 'nSteps' passi e ne calcola la media per passo
std::vector<double> OneBlockMean(RandomWalk* rw, unsigned int throws, unsigned int nSteps);

// Scrive su file i vettori mean e StdDev in due colonne
void outData(std::vector<double> mean, std::vector<double> StdDev, const char *file);

// Calcola progressivamente media e dev.std. su nBlocksMax blocchi,
// ciascuno con nThrows cammini di nSteps passi, salva in 'filename'
void MeanAndDevStd(RandomWalk* rw, unsigned int nBlocksMax,
                   unsigned int nThrows, unsigned int nSteps,
                   const char * filename);

#endif
