// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPROXY_SYSTEM_H
#define COLVARPROXY_SYSTEM_H


/// Methods for accessing the simulation system (PBCs, integrator, etc)
class colvarproxy_system {

public:

  /// Constructor
  colvarproxy_system();

  /// Destructor
  virtual ~colvarproxy_system();

  /// \brief Name of the unit system used internally by Colvars (by default, that of the back-end).
  /// Supported depending on the back-end: real (A, kcal/mol), metal (A, eV), electron (Bohr, Hartree), gromacs (nm, kJ/mol)
  /// Note: calls to back-end PBC functions assume back-end length unit
  /// We use different unit from back-end in VMD bc using PBC functions from colvarproxy base class
  /// Colvars internal units are user specified, because the module exchanges info in unknown
  /// composite dimensions with user input, while it only exchanges quantities of known
  /// dimension with the back-end (length and forces)
  std::string units;

  /// \brief Convert a length from Angstrom to internal
  inline cvm::real angstrom_to_internal(cvm::real l) const
  {
    return l * angstrom_value_;
  }

  /// \brief Convert a length from internal to Angstrom
  inline cvm::real internal_to_angstrom(cvm::real l) const
  {
    return l / angstrom_value_;
  }

  /// Boltzmann constant, with unit the same as energy / K
  inline cvm::real boltzmann() const
  {
    return boltzmann_;
  }

  /// Current target temperature of the simulation (K units)
  inline cvm::real target_temperature() const
  {
    return target_temperature_;
  }

  /// Set the current target temperature of the simulation (K units)
  virtual int set_target_temperature(cvm::real T);

  /// Time step of the simulation (fs units)
  inline double dt() const
  {
    return timestep_;
  }

  /// Set the current integration timestep of the simulation (fs units)
  virtual int set_integration_timestep(cvm::real dt);

  /// \brief Pseudo-random number with Gaussian distribution
  virtual cvm::real rand_gaussian(void);

  /// Pass restraint energy value for current timestep to MD engine
  virtual void add_energy(cvm::real energy);

  /// \brief Tell the proxy whether total forces are needed (they may not
  /// always be available)
  virtual void request_total_force(bool yesno);

  /// Are total forces being used?
  virtual bool total_forces_enabled() const;

  /// Are total forces from the current step available?
  virtual bool total_forces_same_step() const;

  /// Get value of alchemical lambda parameter from back-end (if available)
  virtual int get_alch_lambda(cvm::real* lambda);

  /// Set value of alchemical lambda parameter to be sent to back-end at end of timestep
  void set_alch_lambda(cvm::real lambda);

  /// Send cached value of alchemical lambda parameter to back-end (if available)
  virtual int send_alch_lambda();

  /// Get energy derivative with respect to lambda (if available)
  virtual int get_dE_dlambda(cvm::real* dE_dlambda);

  /// Apply a scalar force on dE_dlambda (back-end distributes it onto atoms)
  virtual int apply_force_dE_dlambda(cvm::real* force);

  /// Get energy second derivative with respect to lambda (if available)
  virtual int get_d2E_dlambda2(cvm::real* d2E_dlambda2);

  /// Force to be applied onto alch. lambda, propagated from biasing forces on dE_dlambda
  cvm::real indirect_lambda_biasing_force;

  /// Get weight factor from accelMD
  virtual cvm::real get_accelMD_factor() const {
    cvm::error("Error: accessing the reweighting factor of accelerated MD  "
               "is not yet implemented in the MD engine.\n",
               COLVARS_NOT_IMPLEMENTED);
    return 1.0;
  }
  virtual bool accelMD_enabled() const {
    return false;
  }

  inline int get_dim() const {
    return n_dim;
  }

  inline cvm::vector1d<cvm::real> get_positions() const {
    return positions;
  }

  cvm::vector1d<cvm::real> colvar_forces;

protected:

  //! Total number of Colvars atoms
  size_t n_dim = 0;

  cvm::vector1d<cvm::real> positions;

  /// Next value of lambda to be sent to back-end
  cvm::real cached_alch_lambda;

  /// Whether lambda has been set and needs to be updated in backend
  bool cached_alch_lambda_changed;

  /// Boltzmann constant in internal Colvars units
  cvm::real boltzmann_;

  /// Most up to date target temperature (K units); default to 0.0 if undefined
  cvm::real target_temperature_;

  /// Current integration timestep (engine units); default to 1.0 if undefined
  double timestep_;

  /// \brief Value of 1 Angstrom in the internal (front-end) Colvars unit for atomic coordinates
  /// * defaults to 0 in the base class; derived proxy classes must set it
  /// * in VMD proxy, can only be changed when no variables are defined
  /// as user-defined values in composite units must be compatible with that system
  cvm::real angstrom_value_;

  /// \brief Value of 1 kcal/mol in the internal Colvars unit for energy
  cvm::real kcal_mol_value_;

  /// Whether the total forces have been requested
  bool total_force_requested;

  /// \brief Type of boundary conditions
  ///
  /// Orthogonal and triclinic cells are made available to objects.
  /// For any other conditions (mixed periodicity, triclinic cells in LAMMPS)
  /// minimum-image distances are computed by the host engine regardless.
  enum Boundaries_type {
    boundaries_non_periodic,
    boundaries_pbc_ortho,
    boundaries_pbc_triclinic,
    boundaries_unsupported
  };

  /// Type of boundary conditions
  Boundaries_type boundaries_type;

};

#endif
