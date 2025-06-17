/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cmath>
#include <cstdlib>
#include <string>
#include "system.h"

using namespace std;
using namespace arma;

void System :: step(){
  this->move(); // Perform a MC step on a randomly choosen particle
  _nattempts ++; //update number of attempts performed on the system
  return;
}


void System :: move(){ // Propose a MC move for particle i 
  vec shift(_ndim);       // Will store the proposed translation
  for(int j=0; j<_ndim; j++){
    shift(j) = _rnd.Rannyu(-1.0,1.0) * _delta; // uniform distribution in [-_delta;_delta)
  }
  double xold = _particle.getposition(0, false);
  _particle.translate(shift);  //Call the function Particle::translate
  double xnew = _particle.getposition(0, true);
  //cout << endl << xold << "  a  " << xnew << endl;
  if(this->metro(xold, xnew)){ //Metropolis acceptance evaluation
    _particle.acceptmove();
    _naccepted++;
  } else _particle.moveback(); //If translation is rejected, restore the old configuration
  //cout << endl << xold << "  b  " << xnew << endl;

  
  return;
}

bool System :: metro(double xold, double xnew){ // Metropolis algorithm
  bool decision = false;
  double acceptance;
  double psi_new_squared=pow(exp(-pow(xnew-_mu,2)/double(2*_sigma*_sigma))+exp(-pow(xnew+_mu,2)/double(2*_sigma*_sigma)),2);
  double psi_old_squared=pow(exp(-pow(xold-_mu,2)/double(2*_sigma*_sigma))+exp(-pow(xold+_mu,2)/double(2*_sigma*_sigma)),2);
  acceptance=psi_new_squared/double(psi_old_squared);
  //cout << endl << "acceptance : " << acceptance << "    sigma: " << _sigma << endl;

  //cout << endl << xnew << "  c  " << xold << endl;
  if(_rnd.Rannyu() < acceptance ){
    decision = true; //Metropolis acceptance step
    _appo_psi=psi_new_squared;
  } else{
    _appo_psi=psi_old_squared;
  }
  return decision;
}

void System :: read_configuration(){
  ifstream cinf;
  cinf.open("../INPUT/CONFIG/config.x");
  if(cinf.is_open()){
    string comment;
    string particle;
    double x;
    int ncoord;
    cinf >> ncoord;
    if (ncoord != 1){
      cerr << "PROBLEM: conflicting number of coordinates in input.dat & config.xyz not match!" << endl;
      exit(EXIT_FAILURE);
    }
    cinf >> comment;
    cinf >> particle >> x; 
    _particle.setposition(0, x);
    _particle.acceptmove(); // _x_old = _x_new
    
  } else cerr << "PROBLEM: Unable to open INPUT file config.x"<< endl;
  cinf.close();
  return;
}

void System :: initialize(){ // Initialize the System object according to the content of the input files in the ../INPUT/ directory

  int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
  ifstream Primes("../INPUT/Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  int seed[4]; // Read the seed of the RNG
  ifstream Seed("../INPUT/seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  _rnd.SetRandom(seed,p1,p2);
  _particle.initialize();

  ofstream couta("../OUTPUT/acceptance.dat"); // Set the heading line in file ../OUTPUT/acceptance.dat
  couta << "#   N_BLOCK:  ACCEPTANCE:" << endl;
  couta.close();

  ifstream input("../INPUT/input.dat"); // Start reading ../INPUT/input.dat
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat");
  string property;
  double delta;
  while ( !input.eof() ){
    input >> property;
    if( property == "RESTART" ){
      input >> _restart;
    } else if( property == "DELTA" ){
      input >> delta;
      coutf << "DELTA= " << delta << endl;
      _delta = delta;
    } else if( property == "MU" ){
      input >> _mu;
      coutf << "MU= " << _mu << endl;
    } else if( property == "SIGMA" ){
      input >> _sigma;
      coutf << "SIGMA= " << _sigma << endl;
    } else if( property == "NBLOCKS" ){
      input >> _nblocks;
      coutf << "NBLOCKS= " << _nblocks << endl;
    } else if( property == "NSTEPS" ){
      input >> _nsteps;
      coutf << "NSTEPS= " << _nsteps << endl;
    } else if( property == "EQSTEPS" ){
      input >> _eqsteps;
      coutf << "EQSTEPS= " << _eqsteps << endl;
    } else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << endl;
      break;
    } else cerr << "PROBLEM: unknown input" << endl;
  }
  input.close();
  this->read_configuration();
  coutf << "System initialized!" << endl;
  coutf.close();
  return;
}


void System :: initialize_properties(){ // Initialize data members used for measurement of properties

  string property;
  int index_property = 0;
  _nprop = 0;

  _measure_penergy  = false; //Defining which properties will be measured
  _measure_kenergy  = false;
  _measure_H_av  = false;
  _measure_psi = false;


  ifstream input("../INPUT/properties.dat");
  if (input.is_open()){
    while ( !input.eof() ){
      input >> property;
      if( property == "ENDPROPERTIES" ){
        ofstream coutf;
        coutf.open("../OUTPUT/output.dat",ios::app);
        coutf << "Reading properties completed!" << endl;
        coutf.close();
        break;
      } else if( property == "H_AVERAGE" ){
        _measure_kenergy = true;
        _measure_penergy = true;
        ofstream coutt("../OUTPUT/H_av.dat");
        coutt << "#     BLOCK:   ACTUAL_H:    H_AVE:      ERROR:" << endl;
        coutt.close();
        _nprop++;
        _measure_H_av = true;
        _index_H_av = index_property;
        index_property++;
      } else if( property == "KINETIC_ENERGY"){
        ofstream coutk("../OUTPUT/kinetic_energy.dat");
        coutk << "#     BLOCK:   ACTUAL_KE:    KE_AVE:      ERROR:" << endl;
        coutk.close();
        _nprop++;
        _measure_kenergy = true;
        _index_kenergy = index_property;
        index_property++;
      } else if ( property == "PSI_SQUARED"){
        ofstream coutp("../OUTPUT/psi_squared.dat");
        coutp << "#     BLOCK:  ACTUAL_PSI:     PSI_AVE:      ERROR:" << endl;
        coutp.close();
        input>>_n_bins_psi;
        _nprop += _n_bins_psi;
        _bin_size_psi = 4./_n_bins_psi; // TO BE FIXED IN EXERCISE 4
        _measure_psi = true;
        _index_psi = index_property;
        index_property += _n_bins_psi;
      } else if ( property == "POTENTIAL_ENERGY"){
        ofstream coutp("../OUTPUT/potential_energy.dat");
        coutp << "#     BLOCK:  ACTUAL_PE:     PE_AVE:      ERROR:" << endl;
        coutp.close();
        _nprop++;
        _index_penergy = index_property;
        _measure_penergy = true;
        index_property++;
      }  else cerr << "PROBLEM: unknown property" << endl;
    }
    input.close();
  } else cerr << "PROBLEM: Unable to open properties.dat" << endl;

  // according to the number of properties, resize the vectors _measurement,_average,_block_av,_global_av,_global_av2
  _measurement.resize(_nprop);
  _average.resize(_nprop);
  _block_av.resize(_nprop);
  _global_av.resize(_nprop);
  _global_av2.resize(_nprop);
  _average.zeros();
  _global_av.zeros();
  _global_av2.zeros();
  _nattempts = 0;
  _naccepted = 0;

  return;
}

void System :: finalize(){
  this->write_configuration();
  _rnd.SaveSeed();
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat",ios::app);
  coutf << "Simulation completed!" << endl;
  coutf.close();
  return;
}

// Write final configurations as .xyz files in the directory ../OUTPUT/CONFIG/
void System :: write_configuration(){
  ofstream coutf;
  coutf.open("../OUTPUT/CONFIG/config.x");
  if(coutf.is_open()){
    coutf << "1" << endl; //Ho solo una particella
    coutf << "#Comment!" << endl;
    coutf << "LJ" << "  " 
          << setprecision(17) << _particle.getposition(0,true) << endl; // x
  } else cerr << "PROBLEM: Unable to open config.x" << endl;
  coutf.close();
  coutf.open("../OUTPUT/CONFIG/conf-1.x");
  if(coutf.is_open()){
    coutf << "1" << endl; //Ho solo una particella
    coutf << "#Comment!" << endl;
    coutf << "LJ" << "  "
          << setprecision(17) << _particle.getposition(0,false)<< endl; // x
    
  } else cerr << "PROBLEM: Unable to open conf-1.xyz" << endl;
  coutf.close();
  return;
}

// Write configuration nconf as a .xyz file in directory ../OUTPUT/CONFIG/
void System :: write_XYZ(int nconf){
  ofstream coutf;
  coutf.open("../OUTPUT/CONFIG/config_" + to_string(nconf) + ".x");
  if(coutf.is_open()){
    coutf << "1" << endl; //Ho solo una particella
    coutf << "#Comment!" << endl;
    coutf << "LJ" << "  " << setw(16) << _particle.getposition(0,true) << endl ;         // x
    
  } else cerr << "PROBLEM: Unable to open config.x" << endl;
  coutf.close();
  return;
}

void System :: equilibrate(){
  cout << endl << "Setting accepted moves counter to 0..." << endl;
  _naccepted=0;
  _nattempts=0;
  for(int j=0; j < _eqsteps; j++){ //loop over steps in a block
    this->step();
    if(j%(_eqsteps/100)==0){
      this->autotuning();
    }
  }
  cout << endl << "Acceptance rate in equilibration phase: " << double(_naccepted)/double(_nattempts) << endl;
  cout << endl << "Setting accepted moves counter to 0..." << endl;
  _naccepted=0;
  _nattempts=0;
  return;
}


void System :: block_reset(int blk){ // Reset block accumulators to zero
  ofstream coutf;
  if(blk>0){
    coutf.open("../OUTPUT/output.dat",ios::app);
    coutf << "Block completed: " << blk << endl;
    coutf.close();
  }
  _naccepted = 0;
  _nattempts = 0;
  _block_av.zeros();
  return;
}

void System :: measure(){ // Measure properties
  _measurement.zeros();
  double H_av_temp=0.0; // temporary accumulator for potential energy

  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
 if (_measure_penergy){
    //cout << _measurement.size() << "    " << _index_penergy << endl;
    _measurement[_index_penergy]=pow(_particle.getposition(0, true),4)-5./2*pow(_particle.getposition(0,true),2);
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  if (_measure_kenergy){
    double mu, sigma;
    mu=_mu;
    sigma=_sigma;
    //cout << endl << _particle.getposition(0, true) << "   " << _mu << endl;
    double c1   = pow((_particle.getposition(0, true)-mu)/double(sigma), 2);        // (x−μ)^2/σ^2
    double c2   = (c1 - 1.) / double(sigma*sigma);      // (c1−1)/σ^2
    double appo1 = exp(-c1/2.) * c2;
   double c3=pow((_particle.getposition(0, true)+mu)/double(sigma),2);
    double c4   = (c3 - 1.)/(sigma*sigma);   // come hai fatto per c2
    double appo2= exp(-c3/2.) * c4;
    double e1=exp(-c1/2.);
    double e2=exp(-c3/2.);
    _measurement[_index_kenergy]=-0.5* (appo1+appo2)/(e1+e2);
  }
  /////////////////H AVERAGE////////////////////
  if (_measure_H_av) {
    _measurement[_index_H_av]=_measurement[_index_kenergy]+_measurement[_index_penergy];
  }
    ///////////////////PSI_SQUARED////////////////////
  if (_measure_psi){
    double cum_sum=0; //somma cumulativa
    int index;
    double x = _particle.getposition(0, true);
    index=static_cast<int>((x+2)/_bin_size_psi);
    //cout << " index: " << index << endl;
    if (index<_n_bins_psi){
      _measurement(_index_psi+index)+=1;
    }

  }



  _block_av += _measurement; //Update block accumulators

  return;
}

void System :: averages(int blk){

  ofstream coutf;
  double average, sum_average, sum_ave2;

  _average     = _block_av / double(_nsteps);
  _global_av  += _average;
  _global_av2 += _average % _average; // % -> element-wise multiplication

  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  if (_measure_penergy){
    coutf.open("../OUTPUT/potential_energy.dat",ios::app);
    average  = _average(_index_penergy);
    sum_average = _global_av(_index_penergy);
    sum_ave2 = _global_av2(_index_penergy);
    coutf << setw(12) << blk 
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  if (_measure_kenergy){
    coutf.open("../OUTPUT/kinetic_energy.dat",ios::app);
    average  = _average(_index_kenergy);
    sum_average = _global_av(_index_kenergy);
    sum_ave2 = _global_av2(_index_kenergy);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // H average //////////////////////////////////////////////////////////////
  if (_measure_H_av){
    coutf.open("../OUTPUT/H_av.dat",ios::app);
    average  = _average(_index_H_av);
    sum_average = _global_av(_index_H_av);
    sum_ave2 = _global_av2(_index_H_av);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  //////////////////////PSI_SQUARED//////////////////
  if (_measure_psi){
    coutf.open("../OUTPUT/psi_squared.dat",ios::app);
    double norm=1./(_bin_size_psi); 
    unsigned int counter =0; 
    for (unsigned int i=0; i<_n_bins_psi; i++){
      double psi_center=_bin_size_psi*i+0.5*_bin_size_psi-2;
      average  = _average(_index_psi+i);
      sum_average = _global_av(_index_psi+i);
      sum_ave2 = _global_av2(_index_psi+i);

      coutf << setw(12) << blk
          << setw(12) << psi_center
          << setw(12) << sum_average/double(blk)*norm
          << setw(12) << this->error(sum_average, sum_ave2, blk)*norm << endl;
    }

    coutf.close();
  }

  // ACCEPTANCE ////////////////////////////////////////////////////////////////
  double fraction;
  coutf.open("../OUTPUT/acceptance.dat",ios::app);
  if(_nattempts > 0) fraction = double(_naccepted)/double(_nattempts);
  else fraction = 0.0; 
  coutf << setw(12) << blk << setw(12) << fraction << endl;
  coutf.close();
  
  return;
}

double System :: error(double acc, double acc2, int blk){
  if(blk <= 1) return 0.0;
  else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
}

int System :: get_nbl(){
  return _nblocks;
}

int System :: get_nsteps(){
  return _nsteps;
}

void System :: autotuning(){
  if(_nattempts > 0){
    double fraction = double(_naccepted)/double(_nattempts);
    ofstream coutf;
    if (fraction<0.4){
      if (fraction<0.01){
        _delta*=0.5;
      coutf.open("../OUTPUT/acceptance.dat",ios::app);
      coutf << endl << "Previous AR: " << fraction << "    New _delta value: " << _delta << endl;
      } else{
      _delta*=0.97;
      coutf.open("../OUTPUT/acceptance.dat",ios::app);
      coutf << endl << "Previous AR: " << fraction << "    New _delta value: " << _delta << endl;
      }
    }
    if(fraction>0.6){
      if (fraction>0.99){
        _delta*=2.;
        coutf.open("../OUTPUT/acceptance.dat",ios::app);
        coutf << endl << "Previous AR: " << fraction << "    New _delta value: " << _delta << endl;
      }else{
        _delta*=1.03;
        coutf.open("../OUTPUT/acceptance.dat",ios::app);
        coutf << endl << "Previous AR: " << fraction << "    New _delta value: " << _delta << endl;
      }
    }
    coutf.close();
  _naccepted = 0;
  _nattempts = 0;

  } 
  
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
