// functions.h

#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <vector>
using namespace std;

// Simulazione S(t) esatta
double St_function(Random* rnd, double S0, double mu, double sigma, double t);

// Simulazione passo Δt di S_t
double St_function_step(Random* rnd, double S_t_i, double mu, double sigma, double DeltaT);

// Iterazione in nStep da t0 a T
double Iterate_S_t_i(Random* rnd, double nStep, double t0, double T,
                     double mu, double sigma, double S0);

// Call e Put payoff
double OneBlockMean(Random* rnd, unsigned int throws,
                    double S0, double mu, double sigma, double t, double K);
double Call(double mu, double T, double K, double ST);
double Put(double mu, double T, double K, double ST);

// Media e dev.std. di call/put (metodo europeo)
void CallAndPut(Random* rnd, unsigned int nBlocksMax, unsigned int nThrows,
                double S0, double mu, double sigma, double t, double K);

// Versione con discretizzazione in nStep
vector<double> OneBlockMean_iter(Random* rnd, unsigned int throws,
                                 double S0, double mu, double sigma,
                                 double t0, double T, double K, double nStep);
void CallAndPut_iter(Random* rnd, unsigned int nBlocksMax, unsigned int nThrows,
                     double S0, double mu, double sigma,
                     double t0, double T, double K, unsigned int nStep);

// Utility per output
void outData(std::vector<double> mean, std::vector<double> StdDev, const char file[]);

#endif
