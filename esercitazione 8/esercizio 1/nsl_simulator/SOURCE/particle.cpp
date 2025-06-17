/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
// particle.cpp
#include <iostream>
#include <math.h>
#include "particle.h"

using namespace std;

void Particle :: initialize(){
   _x.resize(_ndim);     // allocate coordinate vector
   _xold.resize(_ndim);  // allocate old-coordinate vector
   _v.resize(_ndim);     // allocate velocity vector
   return;
}

void Particle :: translate(vec delta){
   for(unsigned int i = 0; i < _ndim; i++){
     _x(i) = _x(i) + delta(i); // shift current position by δ
   }
}

void Particle :: moveback(){
   _x = _xold;           // restore previous coordinates
}

void Particle :: acceptmove(){
   _xold = _x;           // commit current position as old
}

double Particle :: getposition(int dim, bool xnew){
   if(xnew) return _x(dim);     // return current coordinate
   else     return _xold(dim);  // return previous coordinate
}

void Particle :: setposition(int dim, double position){
   _x(dim) = position;   // set current coordinate
   return;
}

void Particle :: setpositold(int dim, double position){
   _xold(dim) = position; // set previous coordinate
   return;
}

double Particle :: getvelocity(int dim){
   return _v(dim);       // return velocity component
}

vec Particle :: getvelocity(){
   return _v;            // return full velocity vector
}

void Particle :: setvelocity(int dim, double velocity){
   _v(dim) = velocity;   // set velocity component
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
