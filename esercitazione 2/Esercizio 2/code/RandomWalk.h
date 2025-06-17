// randomwalk.h

#ifndef __RandomWalk__
#define __RandomWalk__

#include <vector>
#include <iostream>
using namespace std;

// Controlla se 'number' è in [low, high)
bool CheckInterval(double low, double high, double number);

class RandomWalk {
public:
    // Costruttore: riceve il generatore di numeri casuali
    RandomWalk(Random* RandomGenerator);
    
    // Distruttore virtuale per consentire delete corretto
    virtual ~RandomWalk() {};
    
    // Metodi
    inline void GetCoord();                          // Stampa le coord attuali
    double GetNorm();                                // Restituisce la distanza dall’origine
    virtual void step();                             // Esegue un passo (reticolare)
    vector<double> GenerateWalk(unsigned int nSteps); // Genera cammino e restituisce norme

protected:
    vector<double> coord;    // coordinate correnti
    Random* rnd;             // puntatore al generatore
};

// Cammino continuo: eredita da RandomWalk
class RandomWalk_Continuum : public RandomWalk {
public:
    // Costruttore della classe derivata
    RandomWalk_Continuum(Random* RandomGenerator) : RandomWalk(RandomGenerator) {}
    
    // Override del passo: isotropo nel continuo
    virtual void step() override;
};

#endif
