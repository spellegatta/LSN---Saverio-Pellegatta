/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <cmath> // For mathematical functions like pow, exp, fabs, sqrt.
#include <cstdlib> // For general utilities like exit, general purpose functions.
#include <string> // For string manipulation.
#include <algorithm> // For functions like min.
#include "system.h" // Include the header file for the System class.

using namespace std; // Use the standard namespace.
using namespace arma; // Use the Armadillo namespace.

// This function performs a single Monte Carlo step for the system.
void System :: step( ){
  this->move(); // Perform a MC step on a randomly choosen particle.
  _nattempts ++; //update number of attempts performed on the system. Increment the counter for total attempts.
  return;
}

// This function proposes a Monte Carlo move for the particle.
void System :: move(){ // Propose a MC move for particle i 
  vec shift(_ndim);       // Will store the proposed translation. Vector to hold the displacement in each dimension.
  // Generate a random shift for each dimension within the range [-_delta, _delta).
  for(int j=0; j<_ndim; j++){
    shift(j) = _rnd.Rannyu(-1.0,1.0) * _delta; // uniform distribution in [-_delta;_delta).
  }
  // Store the current position before applying the shift, so it can be restored if rejected.
  double xold = _particle.getposition(0, false);
  _particle.translate(shift);   //Call the function Particle::translate. Apply the proposed translation to the particle.
  // Get the new position of the particle after the translation.
  double xnew = _particle.getposition(0, true);
  //cout << endl << xold << "   a   " << xnew << endl; // Debugging output, commented out.
  // Evaluate if the proposed move is accepted using the Metropolis criterion.
  if(this->metro(xold, xnew)){ //Metropolis acceptance evaluation
    _particle.acceptmove(); // If the move is accepted, update the particle's internal 'old' position to the new one.
    _naccepted++; // Increment the counter for accepted moves.

  } else{
    _particle.moveback(); //If translation is rejected, restore the old configuration. If rejected, revert the particle to its previous position.
  }
  //cout << endl << xold << "   b   " << xnew << endl; // Debugging output, commented out.

  
  return;
}

// This function implements the Metropolis algorithm to decide whether to accept a proposed move.
bool System :: metro(double xold, double xnew){ // Metropolis algorithm
  bool decision = false; // Initialize the decision to 'false' (reject by default).
  double acceptance; // Variable to store the acceptance probability.
  // Calculate the squared probability density for the new position using the current trial wave function parameters (_mu_new, _sigma_new).
  double psi_new_squared=pow(exp(-pow(xnew-_mu_new,2)/double(2*_sigma_new*_sigma_new))+exp(-pow(xnew+_mu_new,2)/double(2*_sigma_new*_sigma_new)),2);
  // Calculate the squared probability density for the old position using the current trial wave function parameters.
  double psi_old_squared=pow(exp(-pow(xold-_mu_new,2)/double(2*_sigma_new*_sigma_new))+exp(-pow(xold+_mu_new,2)/double(2*_sigma_new*_sigma_new)),2);
  // Calculate the acceptance ratio.
  acceptance=psi_new_squared/double(psi_old_squared);
  //cout << endl << "acceptance : " << acceptance << "     sigma: " << _sigma_new << endl; // Debugging output, commented out.
  //cout << endl << xnew << "   c   " << xold << endl; // Debugging output, commented out.
  // Compare a random number with the acceptance probability.
  if(_rnd.Rannyu() < acceptance ){
    decision = true; //Metropolis acceptance step. If random number is less than acceptance, accept the move.
    _appo_psi=psi_new_squared; // Store the new psi_squared value if the move is accepted.
  } else{
    _appo_psi=psi_old_squared; // Store the old psi_squared value if the move is rejected.
  }
  return decision; // Return the decision (true for accepted, false for rejected).
}

// This function reads the initial configuration of the system from a file.
void System :: read_configuration(){
  ifstream cinf; // Declare an input file stream.
  cinf.open("../INPUT/CONFIG/config.x"); // Open the configuration file.
  if(cinf.is_open()){ // Check if the file was successfully opened.
    string comment; // Variable to store a comment line.
    string particle; // Variable to store the particle type.
    double x; // Variable to store the particle's x-coordinate.
    int ncoord; // Variable to store the number of coordinates (expected to be 1).
    cinf >> ncoord; // Read the number of coordinates.
    // Error check: if number of coordinates doesn't match 1, terminate.
    if (ncoord != 1){
      cerr << "PROBLEM: conflicting number of coordinates in input.dat & config.xyz not match!" << endl;
      exit(EXIT_FAILURE); // Exit the program due to an error.
    }
    cinf >> comment; // Read the comment line.
    cinf >> particle >> x;  // Read the particle type and its initial x-position.
    _particle.setposition(0, x); // Set the particle's initial current position.
    _particle.acceptmove(); // _x_old = _x_new. Set the particle's old position to be the same as the initial current position.
    
  } else cerr << "PROBLEM: Unable to open INPUT file config.x"<< endl; // Error message if the configuration file cannot be opened.
  cinf.close(); // Close the input file stream.
  return;
}

// This function initializes the System object by reading simulation parameters from various input files.
void System :: initialize(){ // Initialize the System object according to the content of the input files in the ../INPUT/ directory

  int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG. Prime numbers for RNG initialization.
  ifstream Primes("../INPUT/Primes"); // Open the file containing prime numbers.
  Primes >> p1 >> p2 ; // Read the two prime numbers.
  Primes.close(); // Close the Primes file.
  int seed[4]; // Read the seed of the RNG. Array to hold the seed for the random number generator.
  ifstream Seed("../INPUT/seed.in"); // Open the file containing the RNG seed.
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3]; // Read the four seed integers.
  _rnd.SetRandom(seed,p1,p2); // Initialize the random number generator with the seed and primes.
  _particle.initialize(); // Initialize the particle object (resize its internal vectors).

  ofstream couta("../OUTPUT/acceptance.dat"); // Set the heading line in file ../OUTPUT/acceptance.dat. Open file for logging acceptance rates.
  couta << "#    N_BLOCK:    ACCEPTANCE:" << endl; // Write header for the acceptance file.
  couta.close(); // Close the acceptance file.

  ofstream coutt("../OUTPUT/SA.dat"); // Open file for logging Simulated Annealing progress.
  coutt << "T:    H:    H_ERR:    MU:    SIGMA:" << endl; // Write header for the SA file.
  coutt.close(); // Close the SA file.

  ifstream input("../INPUT/input.dat"); // Start reading ../INPUT/input.dat. Open the main input parameters file.
  ofstream coutf; // Declare an output file stream.
  coutf.open("../OUTPUT/output.dat"); // Open the main output log file.
  string property; // Variable to store the name of the parameter being read.
  double delta; // Temporary variable for the delta parameter.
  // Loop to read parameters from the input file until the end.
  while ( !input.eof() ){
    input >> property; // Read the next property name.
    // Check the property name and assign the corresponding value to the system's members.
    if( property == "RESTART" ){
      input >> _restart; // Read the restart flag.
    } else if( property == "DELTA" ){
      input >> delta; // Read the value for delta.
      coutf << "DELTA= " << delta << endl; // Log delta value to output file.
      _delta = delta; // Assign to the system's _delta member.
      cout << endl << "Delta=" << _delta << endl; // Output delta to console.
    } else if( property == "MU" ){
      input >> _mu_new; // Read the initial value for mu (new).
      _mu_old=_mu_new; // Set the old mu to be the same as the new mu initially.
      coutf << "MU= " << _mu_new << endl; // Log mu value.
    } else if( property == "SIGMA" ){
      input >> _sigma_new; // Read the initial value for sigma (new).
      _sigma_old=_sigma_new; // Set the old sigma to be the same as the new sigma initially.
      coutf << "SIGMA= " << _sigma_new << endl; // Log sigma value.
    } else if( property == "TEMP0" ){
      input >> _temp; // Read the initial temperature for Simulated Annealing.
      coutf << "TEMP0= " << _temp << endl; // Log initial temperature.
    } else if( property == "NBLOCKS" ){
      input >> _nblocks; // Read the number of blocks.
      coutf << "NBLOCKS= " << _nblocks << endl; // Log number of blocks.
    } else if( property == "NSTEPS" ){
      input >> _nsteps; // Read the number of steps per block.
      coutf << "NSTEPS= " << _nsteps << endl; // Log number of steps.
    } else if( property == "EQSTEPS" ){
      input >> _eqsteps; // Read the number of equilibration steps.
      coutf << "EQSTEPS= " << _eqsteps << endl; // Log equilibration steps.
    } else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << endl; // Indicate that input reading is complete.
      break; // Exit the loop.
    } else cerr << "PROBLEM: unknown input" << endl; // Error message for unrecognized input property.
  }
  input.close(); // Close the input file.
  this->read_configuration(); // Read the initial particle configuration.
  coutf << "System initialized!" << endl; // Log that the system initialization is complete.
  coutf.close(); // Close the main output log file.
  return;
}

// This function calculates the Boltzmann factor for the acceptance probability in Simulated Annealing.
double System :: Boltzmann(){
  // The Boltzmann factor is min(1, exp(-(H_new - H_old) / T)).
  double P=min(1., exp(-double(_H_new-_H_old)/_temp) );
  return P; // Return the calculated probability.
}

// This function sets the random number generator object for the system.
void System :: set_rnd(Random rnd){
  _rnd=rnd; // Assign the provided Random object to the system's _rnd member.
  return;
}

// This function proposes a new set of (mu, sigma) parameters in the context of Simulated Annealing.
void System::SA_proposal(double rescale_T) {
    // 1) Copia completa dello stato corrente. // Create a full copy of the current system state.
    System trial = *this; // Create a 'trial' system by copying the current system.
    trial.initialize_properties(); // Re-initialize properties for the trial system (resets measurement vectors).
    trial.block_reset(0); // Reset block accumulators for the trial system.
    trial.set_rnd(this->_rnd); // Set the random number generator for the trial system to be the same as the current system's.

    // 2) Genera la proposta di annealing sul trial. // Generate the annealing proposal on the trial system.
    double T_new = this->_temp * rescale_T; // Calculate the new temperature by scaling down the current temperature.
    this->_temp = T_new;         // <— abbasso sempre la temperatura. // Update the system's temperature to the new, lower value.
    trial._temp = T_new; // Set the trial system's temperature to the new value.
    trial._mu_old    = trial._mu_new; // Store the current _mu_new as _mu_old for the trial system.
    trial._sigma_old = trial._sigma_new; // Store the current _sigma_new as _sigma_old for the trial system.

    double cs, cm; // Variables for random changes to sigma (cs) and mu (cm).
    do {
      //cs = _rnd.Rannyu(-0.05, 0.05) * T_new; // Previous proposal for cs (commented out).
      cs=_rnd.Rannyu(-0.1, 0.1)*pow(T_new, 0.1); // Propose a change to sigma based on a scaled random number and new temperature.
      //cs=_rnd.Rannyu(-0.01, 0.01); // Another previous proposal for cs (commented out).

    } while (trial._sigma_new + cs < 1e-1); // Ensure the new sigma remains above a minimum threshold.
    trial._sigma_new = this->_sigma_new + cs; // Update the trial system's sigma.
    //cm = _rnd.Rannyu(-0.05, 0.05) * T_new; // Previous proposal for cm (commented out).
    cm=_rnd.Rannyu(-0.1, 0.1)*pow(T_new, 0.1); // Propose a change to mu based on a scaled random number and new temperature.
    //cm=_rnd.Rannyu(-0.001, 0.001); // Another previous proposal for cm (commented out).
    trial._mu_new = this->_mu_new + cm; // Update the trial system's mu.
    // 3) Calcola H_new e H_new_err con i nuovi parametri (su trial). // Calculate H_new and H_new_err with the new parameters (on trial system).
    trial.run_internal_cycle(); // Run a full simulation cycle with the proposed new (mu, sigma) parameters to evaluate the Hamiltonian.

    // 4) Metropolis SA step. // Perform the Metropolis acceptance step for Simulated Annealing.
    double H_old = this->_H_new; // Get the Hamiltonian from the current system.
    double H_new = trial._H_new; // Get the Hamiltonian from the trial system.
    double DH    = H_new - H_old; // Calculate the change in Hamiltonian.
    double P     = min(1.0, exp(-DH / T_new)); // Calculate the acceptance probability using the Boltzmann factor.
    //cout << endl << "exp: " << exp(-DH / T_new) << " T_new " << T_new << " DH " << DH<< " H new "<< H_new << " H_old " << H_old << endl; // Debugging output, commented out.
    bool accept  = (_rnd.Rannyu() < P); // Decide whether to accept the trial state.

    // 5) Scrivi i dati su file SA.dat. // Write data to SA.dat file.
    ofstream fout("../OUTPUT/SA.dat", ios::app); // Open the SA log file in append mode.
    fout // Write the temperature, Hamiltonian, its error, mu, and sigma, choosing between accepted or old values.
      << setw(12) << (accept ? T_new      : this->_temp)
      << setw(12) << (accept ? H_new      : H_old)
      << setw(12) << (accept ? trial._H_new_err : this->_H_new_err)
      << setw(12) << (accept ? trial._mu_new    : this->_mu_new)
      << setw(12) << (accept ? trial._sigma_new : this->_sigma_new)
      << "\n"; // Newline character.
    fout.close(); // Close the SA log file.

    // 6) Se accettato, sostituisci lo stato principale con trial. // If accepted, replace the main state with the trial state.
    if (accept) {
        *this = trial; // Use the assignment operator to update the current system's state with the trial system's.
        this->set_rnd(trial.get_rnd()); // Ensure the random number generator state is also copied from the accepted trial.
    }
    // se rifiutato, il this rimane com’era. // If rejected, the 'this' object remains as it was.
    this->_rnd=trial.get_rnd(); // Regardless of acceptance, synchronize the main system's RNG with the trial's (which might have advanced during the internal cycle).
}

// This function returns the random number generator object.
Random System :: get_rnd(){
  return _rnd; // Return the internal Random object.
}

// This function runs a complete internal simulation cycle (equilibration, blocks, and measurements) to evaluate properties for a given set of parameters.
void System::run_internal_cycle() {

  this->block_reset(0); // Reset accumulators before starting the internal cycle.
  this->equilibrate(false); // Equilibrate the system. Autotuning is disabled here as it's an internal cycle for evaluation.
  // Loop over the simulation blocks.
  for (int i = 0; i < this->get_nbl(); ++i) {
    // Loop over steps within each block.
    for (int j = 0; j < this->get_nsteps(); ++j) {
      this->step(); // Perform a Monte Carlo step.
      this->measure(); // Measure properties.
    }
    this->averages(i + 1); // Calculate and write averages for the completed block.
    this->block_reset(i + 1); // Reset accumulators for the next block.
  }

  return;
}

// This function initializes data members used for measurement of properties and sets up output files.
void System :: initialize_properties(){ // Initialize data members used for measurement of properties

  string property; // Variable to store the property name read from the file.
  int index_property = 0; // Index to assign to each measured property.
  _nprop = 0; // Counter for the total number of properties to be measured.

  _measure_penergy  = false; //Defining which properties will be measured. Flag to indicate if potential energy should be measured.
  _measure_kenergy  = false; // Flag to indicate if kinetic energy should be measured.
  _measure_H_av   = false; // Flag to indicate if average Hamiltonian should be measured.
  _measure_psi = false; // Flag to indicate if psi squared distribution should be measured.

  ifstream input("../INPUT/properties.dat"); // Open the properties definition file.
  if (input.is_open()){ // Check if the file was successfully opened.
    while ( !input.eof() ){ // Loop through the file until the end.
      input >> property; // Read the next property name.
      // Process different properties based on their names.
      if( property == "ENDPROPERTIES" ){
        ofstream coutf; // Declare an output file stream.
        coutf.open("../OUTPUT/output.dat",ios::app); // Open the main output log file in append mode.
        coutf << "Reading properties completed!" << endl; // Log completion of properties reading.
        coutf.close(); // Close the output file.
        break; // Exit the loop.
      } else if( property == "H_AVERAGE" ){
        _measure_kenergy = true; // H_AVERAGE requires kinetic energy.
        _measure_penergy = true; // H_AVERAGE requires potential energy.
        ofstream coutt("../OUTPUT/H_av.dat"); // Open output file for average Hamiltonian.
        coutt << "#     BLOCK:    ACTUAL_H:     H_AVE:       ERROR:" << endl; // Write header to the file.
        coutt.close(); // Close the file.
        _nprop++; // Increment total property count.
        _measure_H_av = true; // Set flag to measure average Hamiltonian.
        _index_H_av = index_property; // Assign index for average Hamiltonian.
        index_property++; // Increment property index for the next property.
      } else if( property == "KINETIC_ENERGY"){
        ofstream coutk("../OUTPUT/kinetic_energy.dat"); // Open output file for kinetic energy.
        coutk << "#     BLOCK:    ACTUAL_KE:      KE_AVE:       ERROR:" << endl; // Write header to the file.
        coutk.close(); // Close the file.
        _nprop++; // Increment total property count.
        _measure_kenergy = true; // Set flag to measure kinetic energy.
        _index_kenergy = index_property; // Assign index for kinetic energy.
        index_property++; // Increment property index.
      } else if ( property == "POTENTIAL_ENERGY"){
        ofstream coutp("../OUTPUT/potential_energy.dat"); // Open output file for potential energy.
        coutp << "#     BLOCK:  ACTUAL_PE:      PE_AVE:       ERROR:" << endl; // Write header to the file.
        coutp.close(); // Close the file.
        _nprop++; // Increment total property count.
        _index_penergy = index_property; // Assign index for potential energy.
        _measure_penergy = true; // Set flag to measure potential energy.
        index_property++; // Increment property index.
      } else if ( property == "PSI_SQUARED"){
        ofstream coutp("../OUTPUT/psi_squared.dat"); // Open output file for psi squared distribution.
        coutp << "#     BLOCK:  ACTUAL_PSI:      PSI_AVE:       ERROR:" << endl; // Write header to the file.
        coutp.close(); // Close the file.
        input>>_n_bins_psi; // Read the number of bins for the psi squared histogram.
        _nprop += _n_bins_psi; // Add the number of bins to the total property count.
        _bin_size_psi = 4./_n_bins_psi; // TO BE FIXED IN EXERCISE 4. Calculate the size of each bin.
        _measure_psi = true; // Set flag to measure psi squared distribution.
        _index_psi = index_property; // Assign starting index for psi squared bins.
        index_property += _n_bins_psi; // Increment property index by the number of bins.
      }   else cerr << "PROBLEM: unknown property" << endl; // Error message for unknown property.
    }
    input.close(); // Close the input file.
  } else cerr << "PROBLEM: Unable to open properties.dat" << endl; // Error message if properties file cannot be opened.

  // according to the number of properties, resize the vectors _measurement,_average,_block_av,_global_av,_global_av2. Resize all measurement and average vectors based on the total number of properties.
  _measurement.resize(_nprop);
  _average.resize(_nprop);
  _block_av.resize(_nprop);
  _global_av.resize(_nprop);
  _global_av2.resize(_nprop);
  // Initialize all measurement and average vectors to zero.
  _average.zeros();
  _global_av.zeros();
  _global_av2.zeros();
  _nattempts = 0; // Reset total attempts counter.
  _naccepted = 0; // Reset accepted moves counter.
  return;
}

// This function performs finalization steps at the end of the simulation.
void System :: finalize(){
  this->write_configuration(); // Write the final configuration of the system to a file.
  _rnd.SaveSeed(); // Save the current state of the random number generator for reproducibility.
  ofstream coutf; // Declare an output file stream.
  coutf.open("../OUTPUT/output.dat",ios::app); // Open the main output log file in append mode.
  coutf << "Simulation completed!" << endl; // Write a message indicating simulation completion.
  coutf.close(); // Close the output file.
  return;
}

// Write final configurations as .xyz files in the directory ../OUTPUT/CONFIG/
void System :: write_configuration(){
  ofstream coutf; // Declare an output file stream.
  // Write the final accepted configuration.
  coutf.open("../OUTPUT/CONFIG/config.x"); // Open the file for the final configuration.
  if(coutf.is_open()){ // Check if the file was successfully opened.
    coutf << "1" << endl; //Ho solo una particella. Write the number of particles (1 in this case).
    coutf << "#Comment!" << endl; // Write a comment line.
    coutf << "LJ" << "   " // Write particle type (e.g., "LJ" for Lennard-Jones).
          << setprecision(17) << _particle.getposition(0,true) << endl; // x. Write the x-coordinate of the particle with high precision.
  } else cerr << "PROBLEM: Unable to open config.x" << endl; // Error message if file cannot be opened.
  coutf.close(); // Close the file.
  // Write the last known 'old' configuration (which is typically the last accepted one before potential rejection).
  coutf.open("../OUTPUT/CONFIG/conf-1.x"); // Open the file for the previous configuration.
  if(coutf.is_open()){ // Check if the file was successfully opened.
    coutf << "1" << endl; //Ho solo una particella. Write the number of particles.
    coutf << "#Comment!" << endl; // Write a comment line.
    coutf << "LJ" << "   " // Write particle type.
          << setprecision(17) << _particle.getposition(0,false)<< endl; // x. Write the old x-coordinate of the particle.
    
  } else cerr << "PROBLEM: Unable to open conf-1.xyz" << endl; // Error message if file cannot be opened.
  coutf.close(); // Close the file.
  return;
}

// Write configuration nconf as a .xyz file in directory ../OUTPUT/CONFIG/
void System :: write_XYZ(int nconf){
  ofstream coutf; // Declare an output file stream.
  // Open a configuration file with a name based on the block number (nconf).
  coutf.open("../OUTPUT/CONFIG/config_" + to_string(nconf) + ".x");
  if(coutf.is_open()){ // Check if the file was successfully opened.
    coutf << "1" << endl; //Ho solo una particella. Write the number of particles.
    coutf << "#Comment!" << endl; // Write a comment line.
    coutf << "LJ" << "   " << setw(16) << _particle.getposition(0,true) << endl ;         // x. Write particle type and its current x-position, formatted.
    
  } else cerr << "PROBLEM: Unable to open config.x" << endl; // Error message if file cannot be opened.
  coutf.close(); // Close the file.
  return;
}

// This function equilibrates the system by performing Monte Carlo steps.
void System :: equilibrate(bool do_tuning){
  cout << endl << "Setting accepted moves counter to 0..." << endl; // Inform user that counters are being reset.
  _naccepted=0; // Reset accepted moves counter.
  _nattempts=0; // Reset total attempts counter.
  // Loop for the specified number of equilibration steps.
  for(int j=0; j < _eqsteps; j++){ //loop over steps in a block
    this->step(); // Perform a Monte Carlo step.
    // Periodically perform autotuning of the displacement step if 'do_tuning' is true.
    if(j%(_eqsteps/100)==0 && do_tuning){
      this->autotuning();
    }
  }
  // Print the acceptance rate during the equilibration phase.
  cout << endl << "Acceptance rate in equilibration phase: " << double(_naccepted)/double(_nattempts) << endl;
  // If autotuning was performed, print the final _delta value.
  if (do_tuning) cout << endl << "\u0394 = " << _delta << endl; // Unicode for Delta symbol.
  cout << endl << "Setting accepted moves counter to 0..." << endl; // Inform user about resetting counters again for the main simulation.
  _naccepted=0; // Reset accepted moves counter.
  _nattempts=0; // Reset total attempts counter.
  return;
}

// This function resets block accumulators for a new simulation block.
void System :: block_reset(int blk){ // Reset block accumulators to zero
  ofstream coutf; // Declare an output file stream.
  // If it's not the very first block (blk > 0), log the completion of the previous block.
  if(blk>0){
    coutf.open("../OUTPUT/output.dat",ios::app); // Open the main output log file in append mode.
    coutf << "Block completed: " << blk << endl; // Write a message indicating which block was completed.
    coutf.close(); // Close the file.
  }
  _naccepted = 0; // Reset the counter for accepted moves for the current block.
  _nattempts = 0; // Reset the counter for total attempts for the current block.
  _block_av.zeros(); // Set all elements of the block average vector to zero.
  
  return;
}

// This function measures various properties of the system at the current Monte Carlo step.
void System :: measure(){ // Measure properties
  _measurement.zeros(); // Initialize the measurement vector to zeros for the current step.
  double H_av_temp=0.0; // temporary accumulator for potential energy. This variable appears to be unused.

  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  // If potential energy measurement is enabled.
  if (_measure_penergy){
    //cout << _measurement.size() << "     " << _index_penergy << endl; // Debugging output, commented out.
    // Calculate the potential energy using the given potential function (x^4 - (5/2)x^2) and the particle's current position.
    _measurement[_index_penergy]=pow(_particle.getposition(0, true),4)-5./2*pow(_particle.getposition(0,true),2);
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  // If kinetic energy measurement is enabled.
  if (_measure_kenergy){
    double mu, sigma; // Local variables for the trial wave function parameters.
    mu=_mu_new; // Use the *new* mu value.
    sigma=_sigma_new; // Use the *new* sigma value.
    //cout << endl << _particle.getposition(0, true) << "   " << _mu_new << endl; // Debugging output, commented out.
    // Calculate intermediate terms for the kinetic energy derived from the second derivative of the trial wave function (Laplacian).
    double c1    = pow((_particle.getposition(0, true)-mu)/double(sigma), 2);        // (x−μ)^2/σ^2
    double c2    = (c1 - 1.) / double(sigma*sigma);       // (c1−1)/σ^2
    double appo1 = exp(-c1/2.) * c2;
    double c3=pow((_particle.getposition(0, true)+mu)/double(sigma),2);
    double c4    = (c3 - 1.)/(sigma*sigma);   // come hai fatto per c2. Similar calculation for the second term involving +mu.
    double appo2= exp(-c3/2.) * c4;
    double e1=exp(-c1/2.);
    double e2=exp(-c3/2.);
    // The kinetic energy expectation value for the trial wavefunction.
    _measurement[_index_kenergy]=-0.5* (appo1+appo2)/(e1+e2);
  }
  /////////////////H AVERAGE////////////////////
  // If average Hamiltonian measurement is enabled.
  if (_measure_H_av) {
    // The average Hamiltonian is the sum of kinetic and potential energies.
    _measurement[_index_H_av]=_measurement[_index_kenergy]+_measurement[_index_penergy];
  }

  ///////////////////PSI_SQUARED////////////////////
  // If psi squared distribution measurement is enabled.
  if (_measure_psi){
    double cum_sum=0; //somma cumulativa. Cumulative sum (appears unused for this specific measurement).
    int index; // Index of the bin for the current particle position.
    double x = _particle.getposition(0, true); // Get the current x-position of the particle.
    // Calculate the bin index based on the particle's position. The range is assumed to be from -2 to +2.
    index=static_cast<int>((x+2)/_bin_size_psi);
    //cout << " index: " << index << endl; // Debugging output, commented out.
    // If the calculated index is within the valid range of bins.
    if (index<_n_bins_psi){
      _measurement(_index_psi+index)+=1; // Increment the count for the corresponding bin.
    }
  }

  _block_av += _measurement; //Update block accumulators. Add the current step's measurements to the block accumulators.

  return;
}

// This function calculates and writes the block averages and their errors to output files.
void System :: averages(int blk){

  ofstream coutf; // Declare an output file stream.
  double average, sum_average, sum_ave2; // Variables for current block average, cumulative sum of averages, and cumulative sum of squared averages.

  _average      = _block_av / double(_nsteps); // Calculate the average of properties for the current block.
  _global_av    += _average; // Add the current block average to the global cumulative average.
  _global_av2 += _average % _average; // % -> element-wise multiplication. Add the squared current block average to the global squared cumulative average.

  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  // If potential energy measurement is enabled.
  if (_measure_penergy){
    coutf.open("../OUTPUT/potential_energy.dat",ios::app); // Open the potential energy output file in append mode.
    average  = _average(_index_penergy); // Get the average potential energy for the current block.
    sum_average = _global_av(_index_penergy); // Get the global sum of potential energy averages.
    sum_ave2 = _global_av2(_index_penergy); // Get the global sum of squared potential energy averages.
    // Write block number, current block average, cumulative average, and error to file.
    coutf << setw(12) << blk 
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close(); // Close the file.
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  // If kinetic energy measurement is enabled.
  if (_measure_kenergy){
    coutf.open("../OUTPUT/kinetic_energy.dat",ios::app); // Open the kinetic energy output file in append mode.
    average  = _average(_index_kenergy); // Get the average kinetic energy for the current block.
    sum_average = _global_av(_index_kenergy); // Get the global sum of kinetic energy averages.
    sum_ave2 = _global_av2(_index_kenergy); // Get the global sum of squared kinetic energy averages.
    // Write block number, current block average, cumulative average, and error to file.
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close(); // Close the file.
  }
  // H average //////////////////////////////////////////////////////////////
  // If average Hamiltonian measurement is enabled.
  if (_measure_H_av){
    coutf.open("../OUTPUT/H_av.dat",ios::app); // Open the average Hamiltonian output file in append mode.
    average  = _average(_index_H_av); // Get the average Hamiltonian for the current block.
    sum_average = _global_av(_index_H_av); // Get the global sum of Hamiltonian averages.
    sum_ave2 = _global_av2(_index_H_av); // Get the global sum of squared Hamiltonian averages.
    // Write block number, current block average, cumulative average, and error to file.
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close(); // Close the file.
    _H_new=sum_average/double(blk); // Store the current global average of the Hamiltonian.
    _H_new_err=this->error(sum_average, sum_ave2, blk); // Store the current error of the Hamiltonian average.
  }
  //////////////////////PSI_SQUARED//////////////////
  // If psi squared distribution measurement is enabled.
  if (_measure_psi){
    coutf.open("../OUTPUT/psi_squared.dat",ios::app); // Open the psi squared distribution output file in append mode.
    double norm=1./(_bin_size_psi); // Normalization factor for the histogram (to get probability density).
    unsigned int counter =0; // Unused counter.
    // Loop through each bin to write its average and error.
    for (unsigned int i=0; i<_n_bins_psi; i++){
      double psi_center=_bin_size_psi*i+0.5*_bin_size_psi-2; // Calculate the center position of the current bin.
      average  = _average(_index_psi+i); // Get the count for the current bin in the block average.
      sum_average = _global_av(_index_psi+i); // Get the global sum of counts for the current bin.
      sum_ave2 = _global_av2(_index_psi+i); // Get the global sum of squared counts for the current bin.

      // Write block number, bin center, normalized cumulative average, and normalized error to file.
      coutf << setw(12) << blk
            << setw(12) << psi_center
            << setw(12) << sum_average/double(blk)*norm
            << setw(12) << this->error(sum_average, sum_ave2, blk)*norm << endl;
    }

    coutf.close(); // Close the file.
  }
  // ACCEPTANCE ////////////////////////////////////////////////////////////////
  double fraction; // Variable to store the acceptance fraction.
  coutf.open("../OUTPUT/acceptance.dat",ios::app); // Open the acceptance rate output file in append mode.
  // Calculate the acceptance fraction (accepted moves / total attempts).
  if(_nattempts > 0) fraction = double(_naccepted)/double(_nattempts);
  else fraction = 0.0; // If no attempts, acceptance is 0.
  coutf << setw(12) << blk << setw(12) << fraction << endl; // Write block number and acceptance fraction.
  coutf.close(); // Close the file.
  
  return;
}

// This function calculates the statistical error using the blocking method.
double System :: error(double acc, double acc2, int blk){
  // If there's only one block or less, the error cannot be reliably calculated, so return 0.0.
  if(blk <= 1) return 0.0;
  // Otherwise, calculate the error using the formula for block averaging (sqrt of variance of block means / number of blocks).
  else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
}

// This function returns the total number of blocks set for the simulation.
int System :: get_nbl(){
  return _nblocks; // Return the value of _nblocks.
}

// This function returns the number of Monte Carlo steps performed within each block.
int System :: get_nsteps(){
  return _nsteps; // Return the value of _nsteps.
}

// This function performs autotuning of the Monte Carlo displacement step (_delta) based on the acceptance rate.
void System :: autotuning(){
  if(_nattempts > 0){ // Only perform autotuning if there have been any attempted moves.
    double fraction = double(_naccepted)/double(_nattempts); // Calculate the current acceptance rate.
    ofstream coutf; // Declare an output file stream.
    // If the acceptance rate is too low (less than 0.4).
    if (fraction<0.4){
      if (fraction<0.01){ // If the acceptance rate is extremely low (less than 0.01).
        _delta*=0.5; // Halve the displacement step size to increase acceptance.
      coutf.open("../OUTPUT/acceptance.dat",ios::app); // Open the acceptance log file in append mode.
      coutf << endl << "Previous AR: " << fraction << "     New _delta value: " << _delta << endl; // Log the change in _delta.
      } else{ // If the acceptance rate is moderately low (between 0.01 and 0.4).
      _delta*=0.97; // Slightly reduce the displacement step size.
      coutf.open("../OUTPUT/acceptance.dat",ios::app); // Open the acceptance log file in append mode.
      coutf << endl << "Previous AR: " << fraction << "     New _delta value: " << _delta << endl; // Log the change in _delta.
      }
    }
    // If the acceptance rate is too high (greater than 0.6).
    if(fraction>0.6){
      if (fraction>0.99){ // If the acceptance rate is extremely high (greater than 0.99).
        _delta*=2.; // Double the displacement step size to decrease acceptance.
        coutf.open("../OUTPUT/acceptance.dat",ios::app); // Open the acceptance log file in append mode.
        coutf << endl << "Previous AR: " << fraction << "     New _delta value: " << _delta << endl; // Log the change in _delta.
      }else{ // If the acceptance rate is moderately high (between 0.6 and 0.99).
        _delta*=1.03; // Slightly increase the displacement step size.
        coutf.open("../OUTPUT/acceptance.dat",ios::app); // Open the acceptance log file in append mode.
        coutf << endl << "Previous AR: " << fraction << "     New _delta value: " << _delta << endl; // Log the change in _delta.
      }
    }
    coutf.close(); // Close the output file.
  _naccepted = 0; // Reset accepted moves counter after autotuning.
  _nattempts = 0; // Reset total attempts counter after autotuning.

  } 
  
}

// This function returns the current temperature of the system.
double System :: getT(){
  return _temp; // Return the value of _temp.
}

// This function sets the temperature of the system.
void System :: setT(double T){
  _temp=T; // Set the _temp member to the provided temperature value.
  return;
}

// Copy constructor for the System class.
System::System(const System &other)
: _ndim(other._ndim) // Initialize _ndim from the 'other' object.
, _temp(other._temp) // Initialize _temp from the 'other' object.
, _mu_new(other._mu_new) // Initialize _mu_new from the 'other' object.
, _sigma_new(other._sigma_new) // Initialize _sigma_new from the 'other' object.
, _H_new(other._H_new) // Initialize _H_new from the 'other' object.
, _H_new_err(other._H_new_err) // Initialize _H_new_err from the 'other' object.
, _H_old(other._H_old) // Initialize _H_old from the 'other' object.
, _H_old_err(other._H_old_err) // Initialize _H_old_err from the 'other' object.
, _nblocks(other._nblocks) // Initialize _nblocks from the 'other' object.
, _nsteps(other._nsteps) // Initialize _nsteps from the 'other' object.
, _eqsteps(other._eqsteps) // Initialize _eqsteps from the 'other' object.
, _block_av(other._block_av) // Initialize _block_av from the 'other' object (deep copy).
, _global_av(other._global_av) // Initialize _global_av from the 'other' object (deep copy).
, _global_av2(other._global_av2) // Initialize _global_av2 from the 'other' object (deep copy).
, _particle(other._particle) // Initialize _particle from the 'other' object (deep copy for Particle object).
, _delta(other._delta) // Initialize _delta from the 'other' object.
, _mu_old(other._mu_old) // Initialize _mu_old from the 'other' object.
, _sigma_old(other._sigma_old) // Initialize _sigma_old from the 'other' object.
{
    // copy any additional members here. // Placeholder for any other members that need to be copied.
}

// 3) Define assignment operator:. // Assignment operator for the System class.
System& System::operator=(const System &other) {
    if (this == &other) return *this; // Self-assignment check: if assigning to itself, do nothing and return.

    // copy each member: // Copy all member variables from 'other' to 'this' object.
    _H_old            = other._H_old;
    _H_old_err        = other._H_old_err;
    _ndim             = other._ndim;
    _temp             = other._temp;
    _mu_new           = other._mu_new;
    _mu_old           = other._mu_old;
    _sigma_old        = other._sigma_old;
    _sigma_new        = other._sigma_new;
    _H_new            = other._H_new;
    _H_new_err        = other._H_new_err;
    _nblocks          = other._nblocks;
    _nsteps           = other._nsteps;
    _eqsteps          = other._eqsteps;
    _block_av         = other._block_av; // Deep copy for Armadillo vectors.
    _global_av        = other._global_av; // Deep copy for Armadillo vectors.
    _global_av2       = other._global_av2; // Deep copy for Armadillo vectors.
    _particle         = other._particle; // Deep copy for Particle object.
    _delta            = other._delta;
    
    // copy any additional members here. // Placeholder for any other members that need to be copied.

    return *this; // Return a reference to the current object.
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
