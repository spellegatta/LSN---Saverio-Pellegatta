/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <iostream>
#include <math.h> // Include for mathematical functions like pow.
#include "particle.h" // Include the header file for the Particle class.

using namespace std;

// This function initializes the size of the position and velocity vectors for the particle.
void Particle :: initialize(){
    _x.resize(_ndim); // Resize the current position vector to the dimensionality of the system.
    _xold.resize(_ndim); // Resize the old position vector to the dimensionality of the system.
    _v.resize(_ndim); // Resize the velocity vector to the dimensionality of the system.
    return;
}

// This function translates the particle's position by a given delta vector.
void Particle :: translate(const vec &delta){
    for(unsigned int i=0; i<_ndim; i++){ // Loop through each dimension.
      _x(i) = _x(i) + delta(i); // Update the current position by adding the corresponding delta component.
    }
}

// This function restores the particle's position to its previous (old) state.
void Particle :: moveback(){
    _x = _xold; // Assign the old position vector to the current position vector.
}

// This function accepts the current position as the new old position.
void Particle :: acceptmove(){
    _xold = _x; // Assign the current position vector to the old position vector.
}

// This function returns the position of the particle in a specified dimension.
// The 'xnew' boolean determines whether to return the current or old position.
double Particle :: getposition(int dim, bool xnew) const{
    if(xnew) return _x(dim); // Return the current position if xnew is true.
    else return _xold(dim); // Return the old position if xnew is false.
}

// This function sets the particle's position in a specified dimension.
void Particle :: setposition(int dim, double position){
    _x(dim) = position; // Set the current position in the given dimension.
    return;
}

// This function sets the particle's old position in a specified dimension.
void Particle :: setpositold(int dim, double position){
    _xold(dim) = position; // Set the old position in the given dimension.
    return;
}

// This function returns the particle's velocity in a specified dimension.
double Particle :: getvelocity(int dim) const{
    return _v(dim); // Return the velocity component in the given dimension.
}

// This function returns the entire velocity vector of the particle.
vec Particle :: getvelocity() const{
    return _v; // Return the full velocity vector.
}

// This function sets the particle's velocity in a specified dimension.
void Particle :: setvelocity(int dim, double velocity){
    _v(dim) = velocity; // Set the velocity component in the given dimension.
    return;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
