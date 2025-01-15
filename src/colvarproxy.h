// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPROXY_H
#define COLVARPROXY_H

#include "colvarmodule.h"
#include "colvartypes.h"
#include "colvarproxy_io.h"
#include "colvarproxy_system.h"

/// \file colvarproxy.h
/// \brief Colvars proxy classes
///
/// This file declares the class for the object responsible for interfacing
/// Colvars with other codes (MD engines, VMD, Python).  The \link colvarproxy
/// \endlink class is a derivative of multiple classes, each devoted to a
/// specific task (e.g. \link colvarproxy_atoms \endlink to access data for
/// individual atoms).
///
/// To interface to a new MD engine, the simplest solution is to derive a new
/// class from \link colvarproxy \endlink.  Currently implemented are: \link
/// colvarproxy_lammps, \endlink, \link colvarproxy_namd, \endlink, \link
/// colvarproxy_vmd \endlink.



#if defined(_OPENMP)
#include <omp.h>
#else
struct omp_lock_t;
#endif

/// \brief Methods for multiple-replica communication
class colvarproxy_replicas {

public:

  /// Constructor
  colvarproxy_replicas();

  /// Destructor
  virtual ~colvarproxy_replicas();

  /// \brief Indicate if multi-replica support is available and active
  virtual int replica_enabled();

  /// \brief Index of this replica
  virtual int replica_index();

  /// \brief Total number of replicas
  virtual int num_replicas();

  /// \brief Synchronize replica with others
  virtual void replica_comm_barrier();

  /// \brief Receive data from other replica
  virtual int replica_comm_recv(char* msg_data, int buf_len, int src_rep);

  /// \brief Send data to other replica
  virtual int replica_comm_send(char* msg_data, int msg_len, int dest_rep);

};



/// Interface between Colvars and MD engine (GROMACS, LAMMPS, NAMD, VMD...)
///
/// This is the base class: each engine is supported by a derived class.
class colvarproxy
  : public colvarproxy_system,
    public colvarproxy_replicas,
    public colvarproxy_io
{

public:

  /// Pointer to the main object
  colvarmodule *colvars;

  /// Constructor
  colvarproxy();

  /// Destructor
  ~colvarproxy() override;

  inline std::string const &engine_name() const
  {
    return engine_name_;
  }

  bool io_available() override;

  /// Request deallocation of the module (currently only implemented by VMD)
  virtual int request_deletion();

  /// Whether deallocation was requested
  inline bool delete_requested() const
  {
    return b_delete_requested;
  }

  /// \brief Reset proxy state, e.g. requested atoms
  virtual int reset();

  /// (Re)initialize the module
  virtual int parse_module_config();

  /// (Re)initialize required member data (called after the module)
  virtual int setup();

  /// Whether the engine allows to fully initialize Colvars immediately
  inline bool engine_ready() const
  {
    return engine_ready_;
  }

  /// Enqueue new configuration text, to be parsed as soon as possible
  void add_config(std::string const &cmd, std::string const &conf);

  /// Update data required by Colvars module (e.g. read atom positions)
  ///
  /// TODO Break up colvarproxy_namd and colvarproxy_lammps function into these
  virtual int update_input();

  /// Update data based on the results of a Colvars call (e.g. send forces)
  virtual int update_output();

  /// Carry out operations needed before next simulation step is run
  int end_of_step();

  /// Print a message to the main log
  virtual void log(std::string const &message);

  /// Print a message to the main log and/or let the host code know about it
  virtual void error(std::string const &message);

  /// Record error message (used by VMD to collect them after a script call)
  void add_error_msg(std::string const &message);

  /// Retrieve accumulated error messages
  std::string const & get_error_msgs();

  /// As the name says
  void clear_error_msgs();

  /// Whether a simulation is running (warn against irrecovarable errors)
  inline bool simulation_running() const
  {
    return b_simulation_running;
  }

  /// Is the current step a repetition of a step just executed?
  /// This is set to true when the step 0 of a new "run" command is being
  /// executed, regardless of whether a state file has been loaded.
  inline bool simulation_continuing() const
  {
    return b_simulation_continuing;
  }

  /// Called at the end of a simulation segment (i.e. "run" command)
  int post_run();

  /// Convert a version string "YYYY-MM-DD" into an integer
  int get_version_from_string(char const *version_string);

  /// Get the version number (higher = more recent)
  int version_number() const
  {
    return version_int;
  }

protected:

  /// Whether the engine allows to fully initialize Colvars immediately
  bool engine_ready_;

  /// Collected error messages
  std::string error_output;

  /// Whether a simulation is running (warn against irrecovarable errors)
  bool b_simulation_running;

  /// Is the current step a repetition of a step just executed?
  /// This is set to true when the step 0 of a new "run" command is being
  /// executed, regardless of whether a state file has been loaded.
  bool b_simulation_continuing;

  /// Whether the entire module should be deallocated by the host engine
  bool b_delete_requested;

  /// Integer representing the version string (allows comparisons)
  int version_int;

  /// Track which features have been acknowledged during the last run
  size_t features_hash;

protected:

  /// Name of the simulation engine that the derived proxy object supports
  std::string engine_name_ = "standalone";

  /// Queue of config strings or files to be fed to the module
  void *config_queue_;

};


#endif
