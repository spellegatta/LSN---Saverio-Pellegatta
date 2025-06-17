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

// This function performs a single Monte Carlo step for the system.
void System :: step(){
  this->move(); // Perform a MC step on a randomly choosen particle
  _nattempts ++; //update number of attempts performed on the system
  return;
}

// This function proposes a Monte Carlo move for a particle.
void System :: move(){ // Propose a MC move for particle i 
  vec shift(_ndim);       // Will store the proposed translation
  // Generate a random shift for each dimension of the particle's position.
  for(int j=0; j<_ndim; j++){
    shift(j) = _rnd.Rannyu(-1.0,1.0) * _delta; // uniform distribution in [-_delta;_delta)
  }
  // Store the current position of the particle before attempting a move.
  double xold = _particle.getposition(0, false);
  _particle.translate(shift);   //Call the function Particle::translate
  // Get the new position of the particle after the translation.
  double xnew = _particle.getposition(0, true);
  //cout << endl << xold << "   a   " << xnew << endl; // Debugging output, commented out.
  // Evaluate if the proposed move is accepted based on the Metropolis criterion.
  if(this->metro(xold, xnew)){ //Metropolis acceptance evaluation
    _particle.acceptmove(); // If accepted, update the particle's stored position.
    _naccepted++; // Increment the counter for accepted moves.
  } else _particle.moveback(); //If translation is rejected, restore the old configuration
  //cout << endl << xold << "   b   " << xnew << endl; // Debugging output, commented out.

  
  return;
}

// This function implements the Metropolis algorithm to decide whether to accept a proposed move.
bool System :: metro(double xold, double xnew){ // Metropolis algorithm
  bool decision = false; // Initialize the decision to reject the move.
  double acceptance; // Variable to store the acceptance probability.
  // Calculate the squared probability density for the new position (psi_new^2).
  double psi_new_squared=pow(exp(-pow(xnew-_mu,2)/double(2*_sigma*_sigma))+exp(-pow(xnew+_mu,2)/double(2*_sigma*_sigma)),2);
  // Calculate the squared probability density for the old position (psi_old^2).
  double psi_old_squared=pow(exp(-pow(xold-_mu,2)/double(2*_sigma*_sigma))+exp(-pow(xold+_mu,2)/double(2*_sigma*_sigma)),2);
  // Calculate the acceptance ratio (min(1, psi_new^2 / psi_old^2)).
  acceptance=psi_new_squared/double(psi_old_squared);
  //cout << endl << "acceptance : " << acceptance << "     sigma: " << _sigma << endl; // Debugging output, commented out.

  //cout << endl << xnew << "   c   " << xold << endl; // Debugging output, commented out.
  // Compare a random number with the acceptance probability to make the decision.
  if(_rnd.Rannyu() < acceptance ) decision = true; //Metropolis acceptance step
  return decision; // Return the decision (true for accepted, false for rejected).
}

// This function reads the initial configuration of the system from a file.
void System :: read_configuration(){
  ifstream cinf; // Declare an input file stream.
  cinf.open("../INPUT/CONFIG/config.x"); // Open the configuration file.
  if(cinf.is_open()){ // Check if the file was successfully opened.
    string comment; // Variable to store a comment line from the file.
    string particle; // Variable to store particle type (e.g., "LJ").
    double x; // Variable to store the particle's x-coordinate.
    int ncoord; // Variable to store the number of coordinates.
    cinf >> ncoord; // Read the number of coordinates.
    // Check if the number of coordinates matches the expected value (1 for this simulation).
    if (ncoord != 1){
      cerr << "PROBLEM: conflicting number of coordinates in input.dat & config.xyz not match!" << endl;
      exit(EXIT_FAILURE); // Exit if there's a mismatch.
    }
    cinf >> comment; // Read the comment line.
    cinf >> particle >> x;  // Read the particle type and its x-coordinate.
    _particle.setposition(0, x); // Set the particle's initial position.
    _particle.acceptmove(); // _x_old = _x_new // Make the current position the accepted old position.
    
  } else cerr << "PROBLEM: Unable to open INPUT file config.x"<< endl; // Error message if file opening fails.
  cinf.close(); // Close the input file.
  return;
}

// This function initializes the System object by reading parameters from input files.
void System :: initialize(){ // Initialize the System object according to the content of the input files in the ../INPUT/ directory

  int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
  ifstream Primes("../INPUT/Primes"); // Open the Primes file.
  Primes >> p1 >> p2 ; // Read the prime numbers.
  Primes.close(); // Close the Primes file.
  int seed[4]; // Read the seed of the RNG
  ifstream Seed("../INPUT/seed.in"); // Open the seed file.
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3]; // Read the seed values.
  _rnd.SetRandom(seed,p1,p2); // Initialize the random number generator.
  _particle.initialize(); // Initialize the particle object.

  ofstream couta("../OUTPUT/acceptance.dat"); // Set the heading line in file ../OUTPUT/acceptance.dat
  couta << "#    N_BLOCK:    ACCEPTANCE:" << endl; // Write the header for the acceptance file.
  couta.close(); // Close the acceptance file.

  ifstream input("../INPUT/input.dat"); // Start reading ../INPUT/input.dat
  ofstream coutf; // Declare an output file stream.
  coutf.open("../OUTPUT/output.dat"); // Open the output file.
  string property; // Variable to store the property name read from input.
  double delta; // Variable to store the delta parameter.
  // Loop through the input file until the end.
  while ( !input.eof() ){
    input >> property; // Read the property name.
    // Process different properties based on their names.
    if( property == "RESTART" ){
      input >> _restart; // Read the restart flag.
    } else if( property == "DELTA" ){
      input >> delta; // Read the delta value.
      coutf << "DELTA= " << delta << endl; // Write delta to output.
      _delta = delta; // Assign delta to the system's delta parameter.
    } else if( property == "MU" ){
      input >> _mu; // Read the mu value.
      coutf << "MU= " << _mu << endl; // Write mu to output.
    } else if( property == "SIGMA" ){
      input >> _sigma; // Read the sigma value.
      coutf << "SIGMA= " << _sigma << endl; // Write sigma to output.
    } else if( property == "NBLOCKS" ){
      input >> _nblocks; // Read the number of blocks.
      coutf << "NBLOCKS= " << _nblocks << endl; // Write number of blocks to output.
    } else if( property == "NSTEPS" ){
      input >> _nsteps; // Read the number of steps per block.
      coutf << "NSTEPS= " << _nsteps << endl; // Write number of steps to output.
    } else if( property == "EQSTEPS" ){
      input >> _eqsteps; // Read the number of equilibration steps.
      coutf << "EQSTEPS= " << _eqsteps << endl; // Write equilibration steps to output.
    } else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << endl; // Indicate completion of input reading.
      break; // Exit the loop.
    } else cerr << "PROBLEM: unknown input" << endl; // Error for unknown input.
  }
  input.close(); // Close the input file.
  this->read_configuration(); // Read the initial configuration of the system.
  coutf << "System initialized!" << endl; // Confirm system initialization.
  coutf.close(); // Close the output file.
  return;
}

// This function initializes the data members used for measuring properties during the simulation.
void System :: initialize_properties(){ // Initialize data members used for measurement of properties

  string property; // Variable to store the property name read from the file.
  int index_property = 0; // Index to assign to each measured property.
  _nprop = 0; // Counter for the total number of properties to be measured.

  _measure_penergy  = false; //Defining which properties will be measured // Flag to measure potential energy.
  _measure_kenergy  = false; // Flag to measure kinetic energy.
  _measure_H_av   = false; // Flag to measure average Hamiltonian.

  ifstream input("../INPUT/properties.dat"); // Open the properties definition file.
  if (input.is_open()){ // Check if the file was successfully opened.
    while ( !input.eof() ){ // Loop through the file until the end.
      input >> property; // Read the property name.
      // Process different properties based on their names.
      if( property == "ENDPROPERTIES" ){
        ofstream coutf; // Declare an output file stream.
        coutf.open("../OUTPUT/output.dat",ios::app); // Open the main output file in append mode.
        coutf << "Reading properties completed!" << endl; // Indicate completion.
        coutf.close(); // Close the output file.
        break; // Exit the loop.
      } else if( property == "H_AVERAGE" ){
        _measure_kenergy = true; // Set flag to measure kinetic energy (needed for H_AVERAGE).
        _measure_penergy = true; // Set flag to measure potential energy (needed for H_AVERAGE).
        ofstream coutt("../OUTPUT/H_av.dat"); // Open the output file for average Hamiltonian.
        coutt << "#     BLOCK:    ACTUAL_H:     H_AVE:       ERROR:" << endl; // Write header.
        coutt.close(); // Close the file.
        _nprop++; // Increment the count of properties.
        _measure_H_av = true; // Set flag to measure average Hamiltonian.
        _index_H_av = index_property; // Assign index for H_AVERAGE.
        index_property++; // Increment property index.
      } else if( property == "KINETIC_ENERGY"){
        ofstream coutk("../OUTPUT/kinetic_energy.dat"); // Open output file for kinetic energy.
        coutk << "#     BLOCK:    ACTUAL_KE:      KE_AVE:       ERROR:" << endl; // Write header.
        coutk.close(); // Close the file.
        _nprop++; // Increment property count.
        _measure_kenergy = true; // Set flag to measure kinetic energy.
        _index_kenergy = index_property; // Assign index for kinetic energy.
        index_property++; // Increment property index.
      } else if ( property == "POTENTIAL_ENERGY"){
        ofstream coutp("../OUTPUT/potential_energy.dat"); // Open output file for potential energy.
        coutp << "#     BLOCK:  ACTUAL_PE:      PE_AVE:       ERROR:" << endl; // Write header.
        coutp.close(); // Close the file.
        _nprop++; // Increment property count.
        _index_penergy = index_property; // Assign index for potential energy.
        _measure_penergy = true; // Set flag to measure potential energy.
        index_property++; // Increment property index.
      }   else cerr << "PROBLEM: unknown property" << endl; // Error for unknown property.
    }
    input.close(); // Close the properties file.
  } else cerr << "PROBLEM: Unable to open properties.dat" << endl; // Error if properties file cannot be opened.

  // Resize vectors used for measurements and averages based on the number of properties.
  _measurement.resize(_nprop);
  _average.resize(_nprop);
  _block_av.resize(_nprop);
  _global_av.resize(_nprop);
  _global_av2.resize(_nprop);
  // Initialize all measurement and average vectors to zero.
  _average.zeros();
  _global_av.zeros();
  _global_av2.zeros();
  _nattempts = 0; // Reset attempt counter.
  _naccepted = 0; // Reset accepted moves counter.
  return;
}

// This function performs finalization steps before the simulation ends.
void System :: finalize(){
  this->write_configuration(); // Write the final configuration of the system.
  _rnd.SaveSeed(); // Save the current state of the random number generator.
  ofstream coutf; // Declare an output file stream.
  coutf.open("../OUTPUT/output.dat",ios::app); // Open the main output file in append mode.
  coutf << "Simulation completed!" << endl; // Write a message indicating simulation completion.
  coutf.close(); // Close the output file.
  return;
}

// Write final configurations as .xyz files in the directory ../OUTPUT/CONFIG/
void System :: write_configuration(){
  ofstream coutf; // Declare an output file stream.
  // Write the final accepted configuration to 'config.x'.
  coutf.open("../OUTPUT/CONFIG/config.x"); // Open the output file.
  if(coutf.is_open()){ // Check if the file was successfully opened.
    coutf << "1" << endl; //Ho solo una particella // Write the number of particles (1 in this case).
    coutf << "#Comment!" << endl; // Write a comment line.
    coutf << "LJ" << "   " 
          << setprecision(17) << _particle.getposition(0,true) << endl; // x // Write particle type and its current x-position.
  } else cerr << "PROBLEM: Unable to open config.x" << endl; // Error if file cannot be opened.
  coutf.close(); // Close the file.
  // Write the last rejected configuration (or the previous accepted one) to 'conf-1.x'.
  coutf.open("../OUTPUT/CONFIG/conf-1.x"); // Open the output file.
  if(coutf.is_open()){ // Check if the file was successfully opened.
    coutf << "1" << endl; //Ho solo una particella // Write the number of particles.
    coutf << "#Comment!" << endl; // Write a comment line.
    coutf << "LJ" << "   "
          << setprecision(17) << _particle.getposition(0,false)<< endl; // x // Write particle type and its old x-position.
    
  } else cerr << "PROBLEM: Unable to open conf-1.xyz" << endl; // Error if file cannot be opened.
  coutf.close(); // Close the file.
  return;
}

// Write configuration nconf as a .xyz file in directory ../OUTPUT/CONFIG/
void System :: write_XYZ(int nconf){
  ofstream coutf; // Declare an output file stream.
  // Open a new configuration file with a unique name based on 'nconf'.
  coutf.open("../OUTPUT/CONFIG/config_" + to_string(nconf) + ".x");
  if(coutf.is_open()){ // Check if the file was successfully opened.
    coutf << "1" << endl; //Ho solo una particella // Write the number of particles.
    coutf << "#Comment!" << endl; // Write a comment line.
    coutf << "LJ" << "   " << setw(16) << _particle.getposition(0,true) << endl ;         // x // Write particle type and its current x-position.
    
  } else cerr << "PROBLEM: Unable to open config.x" << endl; // Error if file cannot be opened.
  coutf.close(); // Close the file.
  return;
}

// This function performs an equilibration phase for the simulation.
void System :: equilibrate(){
  cout << endl << "Setting accepted moves counter to 0..." << endl; // Inform user about resetting counters.
  _naccepted=0; // Reset accepted moves counter.
  _nattempts=0; // Reset total attempts counter.
  // Loop for the specified number of equilibration steps.
  for(int j=0; j < _eqsteps; j++){ //loop over steps in a block
    this->step(); // Perform a Monte Carlo step.
    // Periodically perform autotuning of the simulation parameters.
    if(j%(_eqsteps/100)==0){
      this->autotuning();
    }
  }
  // Print the acceptance rate during the equilibration phase.
  cout << endl << "Acceptance rate in equilibration phase: " << double(_naccepted)/double(_nattempts) << endl;
  cout << endl << "Setting accepted moves counter to 0..." << endl; // Inform user about resetting counters again.
  _naccepted=0; // Reset accepted moves counter for the main simulation.
  _nattempts=0; // Reset total attempts counter for the main simulation.
  return;
}

// This function resets the accumulators for a new simulation block.
void System :: block_reset(int blk){ // Reset block accumulators to zero
  ofstream coutf; // Declare an output file stream.
  // If it's not the first block, write a message indicating the completion of the previous block.
  if(blk>0){
    coutf.open("../OUTPUT/output.dat",ios::app); // Open the main output file in append mode.
    coutf << "Block completed: " << blk << endl; // Write the block completion message.
    coutf.close(); // Close the file.
  }
  _naccepted = 0; // Reset accepted moves counter for the current block.
  _nattempts = 0; // Reset total attempts counter for the current block.
  _block_av.zeros(); // Set all elements of the block average vector to zero.
  return;
}

// This function measures the properties of the system at the current state.
void System :: measure(){ // Measure properties
  _measurement.zeros(); // Initialize the measurement vector to zero.
  double H_av_temp=0.0; // temporary accumulator for potential energy // This variable seems unused here.

  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  // If potential energy measurement is enabled, calculate and store it.
  if (_measure_penergy){
    //cout << _measurement.size() << "     " << _index_penergy << endl; // Debugging output, commented out.
    // Calculate potential energy based on the particle's current position and a given potential function (x^4 - (5/2)x^2).
    _measurement[_index_penergy]=pow(_particle.getposition(0, true),4)-5./2*pow(_particle.getposition(0,true),2);
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  // If kinetic energy measurement is enabled, calculate and store it.
  if (_measure_kenergy){
    double mu, sigma; // Declare local variables for mu and sigma.
    mu=_mu; // Assign the system's mu to local mu.
    sigma=_sigma; // Assign the system's sigma to local sigma.
    //cout << endl << _particle.getposition(0, true) << "   " << _mu << endl; // Debugging output, commented out.
    // Calculate intermediate terms for kinetic energy calculation.
    double c1    = pow((_particle.getposition(0, true)-mu)/double(sigma), 2);        // (x−μ)^2/σ^2
    double c2    = (c1 - 1.) / double(sigma*sigma);       // (c1−1)/σ^2
    double appo1 = exp(-c1/2.) * c2;
    double c3=pow((_particle.getposition(0, true)+mu)/double(sigma),2);
    double c4    = (c3 - 1.)/(sigma*sigma);   // come hai fatto per c2 // Similar calculation for c2.
    double appo2= exp(-c3/2.) * c4;
    double e1=exp(-c1/2.);
    double e2=exp(-c3/2.);
    // Calculate kinetic energy based on the second derivative of the trial wavefunction.
    _measurement[_index_kenergy]=-0.5* (appo1+appo2)/(e1+e2);
  }
  /////////////////H AVERAGE////////////////////
  // If average Hamiltonian measurement is enabled, calculate it as the sum of kinetic and potential energies.
  if (_measure_H_av) {
    _measurement[_index_H_av]=_measurement[_index_kenergy]+_measurement[_index_penergy];
  }

  _block_av += _measurement; //Update block accumulators // Add the current measurements to the block accumulators.

  return;
}

// This function calculates and writes the averages of measured properties for the current block.
void System :: averages(int blk){

  ofstream coutf; // Declare an output file stream.
  double average, sum_average, sum_ave2; // Variables for current block average, cumulative sum, and squared cumulative sum.

  _average      = _block_av / double(_nsteps); // Calculate the average of properties for the current block.
  _global_av    += _average; // Add the current block average to the global cumulative average.
  _global_av2 += _average % _average; // % -> element-wise multiplication // Add the squared current block average to the global squared cumulative average.

  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  // If potential energy measurement is enabled, write its block average and error to file.
  if (_measure_penergy){
    coutf.open("../OUTPUT/potential_energy.dat",ios::app); // Open the potential energy output file in append mode.
    average  = _average(_index_penergy); // Get the current block average for potential energy.
    sum_average = _global_av(_index_penergy); // Get the global cumulative sum for potential energy.
    sum_ave2 = _global_av2(_index_penergy); // Get the global squared cumulative sum for potential energy.
    // Write block number, current average, cumulative average, and error to file.
    coutf << setw(12) << blk 
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close(); // Close the file.
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  // If kinetic energy measurement is enabled, write its block average and error to file.
  if (_measure_kenergy){
    coutf.open("../OUTPUT/kinetic_energy.dat",ios::app); // Open the kinetic energy output file in append mode.
    average  = _average(_index_kenergy); // Get the current block average for kinetic energy.
    sum_average = _global_av(_index_kenergy); // Get the global cumulative sum for kinetic energy.
    sum_ave2 = _global_av2(_index_kenergy); // Get the global squared cumulative sum for kinetic energy.
    // Write block number, current average, cumulative average, and error to file.
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close(); // Close the file.
  }
  // H average //////////////////////////////////////////////////////////////
  // If average Hamiltonian measurement is enabled, write its block average and error to file.
  if (_measure_H_av){
    coutf.open("../OUTPUT/H_av.dat",ios::app); // Open the average Hamiltonian output file in append mode.
    average  = _average(_index_H_av); // Get the current block average for average Hamiltonian.
    sum_average = _global_av(_index_H_av); // Get the global cumulative sum for average Hamiltonian.
    sum_ave2 = _global_av2(_index_H_av); // Get the global squared cumulative sum for average Hamiltonian.
    // Write block number, current average, cumulative average, and error to file.
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close(); // Close the file.
  }
  // ACCEPTANCE ////////////////////////////////////////////////////////////////
  double fraction; // Variable to store the acceptance fraction.
  coutf.open("../OUTPUT/acceptance.dat",ios::app); // Open the acceptance output file in append mode.
  // Calculate the acceptance fraction (accepted moves / total attempts).
  if(_nattempts > 0) fraction = double(_naccepted)/double(_nattempts);
  else fraction = 0.0; // If no attempts, acceptance is 0.
  coutf << setw(12) << blk << setw(12) << fraction << endl; // Write block number and acceptance fraction.
  coutf.close(); // Close the file.
  
  return;
}

// This function calculates the statistical error using blocking method.
double System :: error(double acc, double acc2, int blk){
  // If only one block or less, error is zero.
  if(blk <= 1) return 0.0;
  // Calculate the statistical error using the formula for block averaging.
  else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
}

// This function returns the total number of blocks for the simulation.
int System :: get_nbl(){
  return _nblocks;
}

// This function returns the number of steps per block for the simulation.
int System :: get_nsteps(){
  return _nsteps;
}

// This function performs autotuning of the _delta parameter based on the acceptance rate.
void System :: autotuning(){
  if(_nattempts > 0){ // Ensure there have been attempts to calculate acceptance rate.
    double fraction = double(_naccepted)/double(_nattempts); // Calculate the current acceptance rate.
    ofstream coutf; // Declare an output file stream.
    // Adjust _delta if the acceptance rate is too low.
    if (fraction<0.4){
      if (fraction<0.01){ // If acceptance rate is extremely low, halve _delta.
        _delta*=0.5;
      coutf.open("../OUTPUT/acceptance.dat",ios::app); // Open acceptance file in append mode.
      coutf << endl << "Previous AR: " << fraction << "     New _delta value: " << _delta << endl; // Log the change.
      } else{ // If acceptance rate is moderately low, slightly reduce _delta.
      _delta*=0.97;
      coutf.open("../OUTPUT/acceptance.dat",ios::app); // Open acceptance file in append mode.
      coutf << endl << "Previous AR: " << fraction << "     New _delta value: " << _delta << endl; // Log the change.
      }
    }
    // Adjust _delta if the acceptance rate is too high.
    if(fraction>0.6){
      if (fraction>0.99){ // If acceptance rate is extremely high, double _delta.
        _delta*=2.;
        coutf.open("../OUTPUT/acceptance.dat",ios::app); // Open acceptance file in append mode.
        coutf << endl << "Previous AR: " << fraction << "     New _delta value: " << _delta << endl; // Log the change.
      }else{ // If acceptance rate is moderately high, slightly increase _delta.
        _delta*=1.03;
        coutf.open("../OUTPUT/acceptance.dat",ios::app); // Open acceptance file in append mode.
        coutf << endl << "Previous AR: " << fraction << "     New _delta value: " << _delta << endl; // Log the change.
      }
    }
    coutf.close(); // Close the output file.
  _naccepted = 0; // Reset accepted moves counter after autotuning.
  _nattempts = 0; // Reset total attempts counter after autotuning.

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
