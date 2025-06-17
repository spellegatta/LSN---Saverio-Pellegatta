/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Particle__
#define __Particle__

#include <armadillo>
#include "random.h"

using namespace std;
using namespace arma;

class Particle {

private:
  const int _ndim = 1; // Dimensionality of the system
  vec _x;               // Current position vector
  vec _xold;            // Previous position vector (used in moveback())
  vec _v;               // Velocity vector

public: // Function declarations
  void initialize();                      // Initialize particle properties
  void translate(vec delta);   // Translate the particle within the simulation box
  void moveback();                       // Move particle back to previous position
  void acceptmove();                     // Accept the proposed move and update particle properties
  double getposition(int dim, bool xnew);// Get the position of the particle along a specific dimension
  void   setposition(int dim, double position); // Set the position of the particle along a specific dimension
  void   setpositold(int dim, double position); // Set the previous position of the particle along a specific dimension
  double getvelocity(int dim);           // Get the velocity of the particle along a specific dimension
  vec    getvelocity();                  // Get the velocity vector of the particle
  void   setvelocity(int dim, double velocity); // Set the velocity of the particle along a specific dimension

};

#endif // __Particle__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
