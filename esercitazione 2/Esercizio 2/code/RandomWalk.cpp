// randomwalk.cpp

#include <vector>
#include "random.h"
#include <cmath>
#include "RandomWalk.h"
using namespace std;

// Costruttore: inizializza coord a (0,0,0) e assegna il generatore rnd
RandomWalk::RandomWalk(Random* RandomGenerator){
    vector<double> x(3,0);
    coord = x;                 // coord[0..2] = 0
    rnd = RandomGenerator;     // salvo il puntatore al generatore di numeri casuali
};

// Passo su reticolo: scelgo una delle 3 direzioni con probabilità uniforme
void RandomWalk::step(){
    double x, appo;
    appo = rnd->Rannyu() * 3;  // appo ∈ [0,3): mi serve per selezionare asse
    x = rnd->Rannyu();         // x ∈ [0,1): mi serve per direzione ±
    unsigned int i = 0;
    while (i < 3){
        // se appo ∈ [i, i+1), uso l'asse i
        if (CheckInterval(i, i+1, appo)){
            // metà dei casi decremento, metà incremento
            if (CheckInterval(0, 0.5, x)){
                coord.at(i) += -1;  // passo negativo
            }
            else{
                coord.at(i) += 1;   // passo positivo
            }
            return;  // esco non appena ho fatto il passo
        }
        i++;
    }
}

// Passo continuo: genero direzione isotropa su sfera unitaria
void RandomWalk_Continuum::step(){
    double theta, phi;
    theta = rnd->Rannyu() * M_PI;      // θ ∈ [0,π)
    phi   = rnd->Rannyu() * 2 * M_PI;  // φ ∈ [0,2π)
    // aggiorno ciascuna componente cartesiana
    coord.at(0) += sin(theta) * cos(phi);
    coord.at(1) += sin(theta) * sin(phi);
    coord.at(2) += cos(theta);
}

// Stampa a video le coordinate correnti (debug)
inline void RandomWalk::GetCoord(){
    cout << endl << "Coordinate del cammino: ";
    for (unsigned int i = 0; i <= 2; i++){
        cout << coord.at(i) << "  ";
    }
    cout << endl;
}

// Genera un cammino di nSteps passi, restituisce vettore delle norme ad ogni passo
vector<double> RandomWalk::GenerateWalk(unsigned int nSteps){
    fill(coord.begin(), coord.end(), 0);  // azzero coordinate all’inizio del cammino
    vector<double> norms(nSteps);
    for (unsigned int i = 0; i < nSteps; i++){
        step();                    // effettuo un passo (reticolare o continuo)
        norms.at(i) = GetNorm();   // calcolo e salvo la distanza dall’origine
    }
    return norms;                  // ritorno il vettore di distanze
}

// Calcola la norma euclidea delle coord correnti
double RandomWalk::GetNorm(){
    double CumSum = 0;
    for (unsigned int i = 0; i < coord.size(); i++){
        CumSum += pow(coord.at(i), 2);  // sommo x²+y²+z²
    }
    return sqrt(CumSum);                // √(x²+y²+z²)
}

// Funzione helper: verifica number ∈ [low, high)
bool CheckInterval(double low, double high, double number){
    if (number < high and number >= low){
        return true;
    }
    else{
        return false;
    }
}
