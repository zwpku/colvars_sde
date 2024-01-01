/// -*- c++ -*-
#ifndef COLVARS_COLVARPROXY_SDE_H
#define COLVARS_COLVARPROXY_SDE_H

#include "colvarmodule.h"
#include "colvarproxy.h"
#include <random>

#include "sde_utils.h"

/// \brief Communication between colvars and the SDE backend (implementation of
/// \link colvarproxy \endlink)
class colvarproxy_sde : public colvarproxy {

protected:
  bool first_timestep;

  std::string config_file;
  size_t restart_frequency_s;
  int64_t previous_step;
  double bias_energy;

  // random number generation.
  std::default_random_engine rng;   
  std::normal_distribution<double> *normal_distribution;

public:
  colvarproxy_sde();
  ~colvarproxy_sde();

  // Initialize colvars.
  void init(t_inputrec *sde_inp, int64_t step, const std::string &prefix, const std::string& filename_config);

  // Called each step before evaluating the force provider
  void update_data(int64_t const step);
  /*! \brief
    * Computes forces.
    */

  virtual void calculateForces( std::vector<double> &x, std::vector<double>& bf);

  void add_energy (cvm::real energy);
  void finish();

  // **************** SYSTEM-WIDE PHYSICAL QUANTITIES ****************
  cvm::real rand_gaussian();
  // **************** SIMULATION PARAMETERS ****************
  size_t restart_frequency();
  std::string restart_output_prefix();
  std::string output_prefix();

  // **************** INPUT/OUTPUT ****************
  /// Print a message to the main log
  void log (std::string const &message);
  /// Print a message to the main log and let the rest of the program handle the error
  void error (std::string const &message);
  /// Print a message to the main log and exit with error code
  void fatal_error (std::string const &message);
  /// Print a message to the main log and exit normally
  void exit (std::string const &message);

};

#endif
