/// -*- c++ -*-

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <iostream>

#include "colvarproxy_sde.h"
#include "colvarproxy_sde_version.h"
#include "sde_utils.h"


//************************************************************
// colvarproxy_sde
colvarproxy_sde::colvarproxy_sde() : colvarproxy() {}

// Colvars Initialization
void colvarproxy_sde::init(t_inputrec *ir, int64_t step,
                               const std::string &prefix,
                               const std::string &filename_config) {


  // Initialize colvars.
  first_timestep = true;
  restart_frequency_s = 0;

  angstrom_value_ = 1.0;

  // Retrieve the number of colvar atoms
  n_dim=ir->dim;

  positions.resize(n_dim);
  colvar_forces.resize(n_dim);

  // Get the thermostat temperature.
  // NOTE: Considers only the first temperature coupling group!
  set_target_temperature(ir->ref_t);

  // random number generation.
  // Seed with the mdp parameter ld_seed, the Langevin dynamics seed.
  rng.seed(ir->ld_seed);
  normal_distribution = new std::normal_distribution<double>(0.0, 1.0);

  /// Handle input filenames and prefix/suffix for colvars files.
  ///
  /// filename_config is the colvars configuration file collected from "-colvars" option.
  /// The output prefix will be the prefix of Gromacs log filename.
  /// or "output" otherwise.
  ///
  /// For restart, 'filename_restart' is the colvars input file for restart,
  /// set by the "-cv_restart" option. It will be NULL otherwise.
  ///

  if(!prefix.empty())
  {
    output_prefix_str = prefix;
  }
  else {
    output_prefix_str = "output";
  }

  restart_output_prefix_str = prefix + ".restart";

  // Retrieve masses and charges from input file
  updated_masses_ = updated_charges_ = true;

  // Get timestep 
  set_integration_timestep(ir->delta_t);

  // initiate module: this object will be the communication proxy
  // colvarmodule pointer is only defined on the Master due to the static pointer to colvarproxy.
  colvars = new colvarmodule(this);

  version_int = get_version_from_string(COLVARPROXY_VERSION);

  if (cvm::debug()) {
    log("Initializing the colvars proxy object.\n");
  }

  cvm::log("Using SDE interface, version "+
    cvm::to_str(COLVARPROXY_VERSION)+".\n");

  add_config("configfile", filename_config.c_str());

  colvarproxy::parse_module_config();
  colvars->update_engine_parameters();
  colvars->setup_input();

  colvars->setup_output();

  if (step != 0) {
    cvm::log("Initializing step number to "+cvm::to_str(step)+".\n");
  }

  colvars->it = colvars->it_restart = step;


} // End colvars initialization.


colvarproxy_sde::~colvarproxy_sde()
{}

void colvarproxy_sde::finish()
{
  colvars->write_restart_file(output_prefix_str+".colvars.state");
  colvars->write_output_files();
}

cvm::real colvarproxy_sde::rand_gaussian()
{
  return  (*normal_distribution)(rng);
}

size_t colvarproxy_sde::restart_frequency()
{
  return restart_frequency_s;
}

void colvarproxy_sde::log (std::string const &message)
{
  std::istringstream is(message);
  std::string line;
  while (std::getline(is, line))
    // Gromacs prints messages on the stderr FILE
    fprintf(stderr, "colvars: %s\n", line.c_str());
}

void colvarproxy_sde::error (std::string const &message)
{
  fatal_error (message);
}

void colvarproxy_sde::fatal_error (std::string const &message)
{
  log(message);
  if (!cvm::debug())
    log("If this error message is unclear, "
	"try recompiling with -DCOLVARS_DEBUG.\n");
  cvm::error("Error in collective variables module.\n");
}

void colvarproxy_sde::exit (std::string const &message)
{
  cvm::error("exit: " + message + "\n");
}


void colvarproxy_sde::update_data(int64_t const step)
{

  if(cvm::debug()) {
    cvm::log(cvm::line_marker);
    cvm::log("colvarproxy_sde, step no. "+cvm::to_str(colvars->it)+"\n"+
	    "Updating internal data.\n");
  }

  // step update on master only due to the call of colvars pointer.
  if (first_timestep) {
    first_timestep = false;
  } else {
    if ( step - previous_step == 1 )
      colvars->it++;
    // Other cases?
  }

  previous_step = step;
}


void colvarproxy_sde::calculateForces( std::vector<double> &x, std::vector<double>& bf)
{

  // Zero the forces on the atoms, so that they can be accumulated by the colvars.
  for (size_t i = 0; i < n_dim; i++) {
    colvar_forces[i] = 0.0;
  }

  // Get the positions 
  for (size_t i = 0; i < n_dim; i++) {
    positions[i] = x[i];
  }

  bias_energy = 0.0;
  // Call the collective variable module to fill atoms_new_colvar_forces
  if (colvars->calc() != COLVARS_OK) {
    cvm::error("Error calling colvars->calc()\n");
  }

  // Pass the applied forces to backend
  for (int i = 0; i < n_dim; i++)
    bf[i] = colvar_forces[i];

  return;
}


// Pass restraint energy value for current timestep to MD engine
void colvarproxy_sde::add_energy (cvm::real energy)
{
  bias_energy += energy;
}

