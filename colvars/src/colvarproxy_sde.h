/// -*- c++ -*-
#ifndef COLVARS_COLVARPROXY_SDE_H
#define COLVARS_COLVARPROXY_SDE_H

#include "colvarmodule.h"
#include "colvaratoms.h"
#include "colvarproxy.h"
#include <random>

#include "sde_utils.h"

typedef double rvec[3];

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

  // Node-local bookkepping data
  //! Total number of Colvars atoms
  int        n_colvars_atoms = 0;
  //! Global indices of the Colvars atoms.
  int       *ind = nullptr;
  //! Unwrapped positions for all Colvars atoms, communicated to all nodes.
  rvec      *x_colvars_unwrapped = nullptr;
  //! Old positions for all Colvars atoms on master.
  rvec      *xa_old_whole = nullptr;
  //! Bias forces on all Colvars atoms
  rvec      *f_colvars = nullptr;
public:
  friend class cvm::atom;
  colvarproxy_sde();
  ~colvarproxy_sde();

  // Initialize colvars.
  void init(t_inputrec *sde_inp, int64_t step, const std::string &prefix, const std::string& filename_config);

  // Called each step before evaluating the force provider
  void update_data(int64_t const step);
  /*! \brief
    * Computes forces.
    *
    * \param[in]    forceProviderInput    struct that collects input data for the force providers
    * \param[in,out] forceProviderOutput   struct that collects output data of the force providers
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
  // **************** PERIODIC BOUNDARY CONDITIONS ****************
  cvm::rvector position_distance (cvm::atom_pos const &pos1,
                                  cvm::atom_pos const &pos2) const;

  // **************** INPUT/OUTPUT ****************
  /// Print a message to the main log
  void log (std::string const &message);
  /// Print a message to the main log and let the rest of the program handle the error
  void error (std::string const &message);
  /// Print a message to the main log and exit with error code
  void fatal_error (std::string const &message);
  /// Print a message to the main log and exit normally
  void exit (std::string const &message);
  /// Request to set the units used internally by Colvars
  int set_unit_system(std::string const &units_in, bool colvars_defined);
  /// Read atom identifiers from a file \param filename name of
  /// the file (usually a PDB) \param atoms array to which atoms read
  /// from "filename" will be appended \param pdb_field (optiona) if
  /// "filename" is a PDB file, use this field to determine which are
  /// the atoms to be set
  int load_atoms (char const *filename,
                           std::vector<cvm::atom> &atoms,
                           std::string const &pdb_field,
                           double const pdb_field_value = 0.0);
  /// Load the coordinates for a group of atoms from a file
  /// (usually a PDB); if "pos" is already allocated, the number of its
  /// elements must match the number of atoms in "filename"
  int load_coords (char const *filename,
                            std::vector<cvm::atom_pos> &pos,
                            const std::vector<int> &indices,
                            std::string const &pdb_field,
                            double const pdb_field_value = 0.0);

  int init_atom(int atom_number);

  int check_atom_id(int atom_number);
  void update_atom_properties(int index);
};

#endif
