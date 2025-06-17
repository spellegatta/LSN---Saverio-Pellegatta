/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
// In particle.h

#ifndef PARTICLE_H
#define PARTICLE_H

#include <armadillo>
using arma::vec;

class Particle {
public:
  // Particle interface
  void initialize();
  void translate(const vec &delta);
  void moveback();
  void acceptmove();
  double getposition(int dim, bool xnew) const;
  void   setposition(int dim, double position);
  void   setpositold(int dim, double position);
  double getvelocity(int dim) const;
  vec    getvelocity() const;
  void   setvelocity(int dim, double velocity);

  Particle()
    : _x(_ndim), _xold(_ndim), _v(_ndim)
  {}

  // Copy constructor (must initialize const _ndim)
  Particle(const Particle &other)
    : _ndim(other._ndim)
    , _x(other._x)
    , _xold(other._xold)
    , _v(other._v)
  {}

  // Assignment operator (cannot reassign const _ndim)
  Particle& operator=(const Particle &other) {
    if (this != &other) {
      _x    = other._x;
      _xold = other._xold;
      _v    = other._v;
    }
    return *this;
  }


private:
  const int _ndim = 1; // Dimensionality of the system
  vec _x;               // Current position vector
  vec _xold;            // Previous position vector (used in moveback())
  vec _v;               // Velocity vector
};

#endif // PARTICLE_H

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
