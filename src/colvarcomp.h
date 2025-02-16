// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARCOMP_H
#define COLVARCOMP_H

#include "potentials.h"

// Declaration of colvar::cvc base class and derived ones.
//
// Future cvc's could be declared on additional header files.
// After the declaration of a new derived class, its metric
// functions must be reimplemented as well.
// If the new cvc has no symmetry or periodicity,
// this can be done straightforwardly by using the macro:
// simple_scalar_dist_functions (derived_class)

#include <memory>

#include "colvarmodule.h"
#include "colvar.h"

#ifdef COLVARS_TORCH
#include <torch/torch.h>
#include <torch/script.h>
#endif


/// \brief Colvar component (base class for collective variables)
///
/// A \link colvar::cvc \endlink object (or an object of a
/// cvc-derived class) implements the calculation of a collective
/// variable, its gradients and any other related physical quantities
/// that depend on microscopic degrees of freedom.
///
/// No restriction is set to what kind of calculation a \link colvar::cvc \endlink
/// object performs (usually an analytical function of atomic coordinates).
/// The only constraints are that: \par
///
/// - The value is calculated by the \link calc_value() \endlink
///   method, and is an object of \link colvarvalue \endlink class.  This
///   provides a transparent way to treat scalar and non-scalar variables
///   alike, and allows an automatic selection of the applicable algorithms.
///
/// - The object provides an implementation \link apply_force() \endlink to
///   apply forces to atoms.  Typically, one or more \link colvarmodule::atom_group
///   \endlink objects are used, but this is not a requirement for as long as
///   the \link colvar::cvc \endlink object communicates with the simulation program.
///
/// <b> If you wish to implement a new collective variable component, you
/// should write your own class by inheriting directly from \link
/// colvar::cvc \endlink, or one of its derived classes (for instance,
/// \link colvar::distance \endlink is frequently used, because it provides
/// useful data and function members for any colvar based on two
/// atom groups).</b>
///
/// The steps are: \par
/// 1. Declare the new class as a derivative of \link colvar::cvc \endlink
///    in the file \link colvarcomp.h \endlink
/// 2. Implement the new class in a file named colvarcomp_<something>.cpp
/// 3. Declare the name of the new class inside the \link colvar \endlink class
///    in \link colvar.h \endlink (see "list of available components")
/// 4. Add a call for the new class in colvar::init_components()
////   (file: colvar.cpp)
///

class colvar::cvc
  : public colvarparse, public colvardeps
{
public:

  /// \brief The name of the object (helps to identify this
  /// cvc instance when debugging)
  std::string name;

  /// \brief Description of the type of collective variable
  ///
  /// Normally this string is set by the parent \link colvar \endlink
  /// object within its constructor, when all \link colvar::cvc \endlink
  /// objects are initialized; therefore the main "config string"
  /// constructor does not need to define it.  If a \link colvar::cvc
  /// \endlink is initialized and/or a different constructor is used,
  /// this variable definition should be set within the constructor.
  std::string function_type;

  /// Keyword used in the input to denote this CVC
  std::string config_key;

  /// \brief Coefficient in the polynomial combination (default: 1.0)
  cvm::real sup_coeff;
  /// \brief Exponent in the polynomial combination (default: 1)
  int       sup_np;

  /// \brief Period of the values of this CVC (default: 0.0, non periodic)
  cvm::real period;

  /// \brief If the component is periodic, wrap around this value (default: 0.0)
  cvm::real wrap_center;

  /// \brief Constructor
  ///
  /// Calls the init() function of the class
  cvc(std::string const &conf);

  /// Current initialization state; TODO remove this when using init() after default constructor
  int init_code = COLVARS_OK;

  /// Set the value of \link function_type \endlink and its dependencies
  int set_function_type(std::string const &type);

  /// An init function should be defined for every class inheriting from cvc
  /// \param conf Contents of the configuration file pertaining to this \link
  /// cvc \endlink
  virtual int init(std::string const &conf);

  /// \brief Initialize dependency tree
  virtual int init_dependencies();

  /// \brief Parse options pertaining to total force calculation
  virtual int init_total_force_params(std::string const &conf);

  /// \brief After construction, set data related to dependency handling
  int setup();

  /// \brief Default constructor (used when \link colvar::cvc \endlink
  /// objects are declared within other ones)
  cvc();

  /// Destructor
  virtual ~cvc();

  /// \brief Implementation of the feature list for colvar
  static std::vector<feature *> cvc_features;

  /// \brief Implementation of the feature list accessor for colvar
  virtual const std::vector<feature *> &features() const
  {
    return cvc_features;
  }
  virtual std::vector<feature *> &modify_features()
  {
    return cvc_features;
  }
  static void delete_features() {
    for (size_t i=0; i < cvc_features.size(); i++) {
      delete cvc_features[i];
    }
    cvc_features.clear();
  }

  /// \brief Obtain data needed for the calculation for the backend
  virtual void read_data();

  /// \brief Calculate the variable
  virtual void calc_value() = 0;

  /// \brief Calculate the atomic gradients, to be reused later in
  /// order to apply forces
  virtual void calc_gradients() {}

  /// \brief Calculate finite-difference gradients alongside the analytical ones, for each Cartesian component
  virtual void debug_gradients();

  /// \brief Calculate the total force from the system using the
  /// inverse atomic gradients
  virtual void calc_force_invgrads();

  /// \brief Calculate the divergence of the inverse atomic gradients
  virtual void calc_Jacobian_derivative();


  /// \brief Return the previously calculated value
  colvarvalue const & value() const;

  /// \brief Return the previously calculated total force
  colvarvalue const & total_force() const;

  /// \brief Return the previously calculated divergence of the
  /// inverse atomic gradients
  colvarvalue const & Jacobian_derivative() const;

  /// \brief Apply the collective variable force, by communicating the
  /// atomic forces to the simulation program (\b Note: the \link ft
  /// \endlink member is not altered by this function)
  ///
  /// Note: multiple calls to this function within the same simulation
  /// step will add the forces altogether \param cvforce The
  /// collective variable force, usually coming from the biases and
  /// eventually manipulated by the parent \link colvar \endlink
  /// object
  virtual void apply_force(colvarvalue const &cvforce) = 0;

  /// \brief Square distance between x1 and x2 (can be redefined to
  /// transparently implement constraints, symmetries and
  /// periodicities)
  ///
  /// colvar::cvc::dist2() and the related functions are
  /// declared as "const" functions, but not "static", because
  /// additional parameters defining the metrics (e.g. the
  /// periodicity) may be specific to each colvar::cvc object.
  ///
  /// If symmetries or periodicities are present, the
  /// colvar::cvc::dist2() should be redefined to return the
  /// "closest distance" value and colvar::cvc::dist2_lgrad(),
  /// colvar::cvc::dist2_rgrad() to return its gradients.
  ///
  /// If constraints are present (and not already implemented by any
  /// of the \link colvarvalue \endlink types), the
  /// colvar::cvc::dist2_lgrad() and
  /// colvar::cvc::dist2_rgrad() functions should be redefined
  /// to provide a gradient which is compatible with the constraint,
  /// i.e. already deprived of its component normal to the constraint
  /// hypersurface.
  ///
  /// Finally, another useful application, if you are performing very
  /// many operations with these functions, could be to override the
  /// \link colvarvalue \endlink member functions and access directly
  /// its member data.  For instance: to define dist2(x1,x2) as
  /// (x2.real_value-x1.real_value)*(x2.real_value-x1.real_value) in
  /// case of a scalar \link colvarvalue \endlink type.
  virtual cvm::real dist2(colvarvalue const &x1,
                          colvarvalue const &x2) const;

  /// \brief Gradient(with respect to x1) of the square distance (can
  /// be redefined to transparently implement constraints, symmetries
  /// and periodicities)
  virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;

  /// \brief Gradient(with respect to x2) of the square distance (can
  /// be redefined to transparently implement constraints, symmetries
  /// and periodicities)
  virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
                                  colvarvalue const &x2) const;

  /// \brief Wrap value (for periodic/symmetric cvcs)
  virtual void wrap(colvarvalue &x_unwrapped) const;

  /// Pointer to the gradient of parameter param_name
  virtual colvarvalue const *get_param_grad(std::string const &param_name);

  /// Set the named parameter to the given value
  virtual int set_param(std::string const &param_name, void const *new_value);

  /// \brief Whether or not this CVC will be computed in parallel whenever possible
  bool b_try_scalable;

  /// Forcibly set value of CVC - useful for driving an external coordinate,
  /// eg. lambda dynamics
  inline void set_value(colvarvalue const &new_value) {
    x = new_value;
  }

protected:

  /// Update the description string based on name and type
  int update_description();

  cvm::vector1d<cvm::real> pos;

  /// Record the type of this class as well as those it is derived from
  std::vector<std::string> function_types;

  /// \brief Cached value
  colvarvalue x;

  /// \brief Value at the previous step
  colvarvalue x_old;

  /// \brief Calculated total force (\b Note: this is calculated from
  /// the total atomic forces read from the program, subtracting fromt
  /// the "internal" forces of the system the "external" forces from
  /// the colvar biases)
  colvarvalue ft;

  /// \brief Calculated Jacobian derivative (divergence of the inverse
  /// gradients): serves to calculate the phase space correction
  colvarvalue jd;

  /// \brief Set two scalar boundaries (convenience function)
  void init_scalar_boundaries(cvm::real lb, cvm::real ub);

  /// \brief Location of the lower boundary (not defined by user choice)
  colvarvalue lower_boundary;

  /// \brief Location of the upper boundary (not defined by user choice)
  colvarvalue upper_boundary;

  /// \brief CVC-specific default colvar width
  cvm::real width;
};


inline colvarvalue const & colvar::cvc::value() const
{
  return x;
}


inline colvarvalue const & colvar::cvc::total_force() const
{
  return ft;
}


inline colvarvalue const & colvar::cvc::Jacobian_derivative() const
{
  return jd;
}

class colvar::componentDisabled
  : public colvar::cvc
{
public:
    componentDisabled(std::string const & /* conf */) {
        cvm::error("Error: this component is not enabled in the current build; please see https://colvars.github.io/README-c++11.html");
    }
    virtual ~componentDisabled() {}
    virtual void calc_value() {}
    virtual void calc_gradients() {}
    virtual void apply_force(colvarvalue const & /* force */) {}
};

class colvar::coordinate
  : public colvar::cvc
{
protected:
    // the index of coordinate component
    size_t index;
    cvm::vector1d<cvm::real> grad;

public:
    coordinate(std::string const &conf);
    virtual ~coordinate();
    virtual void calc_value();
    virtual void calc_gradients();
    virtual void apply_force(colvarvalue const &force);
};

class colvar::empiricalcv
  : public colvar::cvc
{
protected:
    cvm::vector1d<cvm::real> grad;
    generic_cv * empirical_cv;
public:
    empiricalcv(std::string const &conf);
    virtual ~empiricalcv();
    virtual void calc_value();
    virtual void calc_gradients();
    virtual void apply_force(colvarvalue const &force);
};


#ifdef COLVARS_TORCH 
// only when LibTorch is available
class colvar::torchANN
  : public colvar::cvc
{
protected:
    torch::jit::script::Module nn;
    /// the index of nn output component
    size_t m_output_index;
    bool use_double_input;
    //bool use_gpu;
    // 1d tensor, concatenation of values of sub-cvcs
    torch::Tensor input_tensor;
    torch::Tensor nn_outputs;

    cvm::vector1d<cvm::real> grad;

public:
    torchANN(std::string const &conf);
    virtual ~torchANN();
    virtual void calc_value();
    virtual void calc_gradients();
    virtual void apply_force(colvarvalue const &force);

    /// Redefined to handle periodicity
    virtual cvm::real dist2(colvarvalue const &x1,
			    colvarvalue const &x2) const;
    /// Redefined to handle periodicity
    virtual colvarvalue dist2_lgrad(colvarvalue const &x1,
				    colvarvalue const &x2) const;
    /// Redefined to handle periodicity
    virtual colvarvalue dist2_rgrad(colvarvalue const &x1,
				    colvarvalue const &x2) const;
    /// Redefined to handle periodicity
    virtual void wrap(colvarvalue &x_unwrapped) const;
};
#else
class colvar::torchANN
  : public colvar::componentDisabled
{
public:
    torchANN(std::string const &conf) : componentDisabled(conf) {}
};
#endif // COLVARS_TORCH checking

// metrics functions for cvc implementations

// simple definitions of the distance functions; these are useful only
// for optimization (the type check performed in the default
// colvarcomp functions is skipped)

// definitions assuming the scalar type

#define simple_scalar_dist_functions(TYPE)                              \
                                                                        \
                                                                        \
  cvm::real colvar::TYPE::dist2(colvarvalue const &x1,                  \
                                colvarvalue const &x2) const            \
  {                                                                     \
    return (x1.real_value - x2.real_value)*(x1.real_value - x2.real_value); \
  }                                                                     \
                                                                        \
                                                                        \
  colvarvalue colvar::TYPE::dist2_lgrad(colvarvalue const &x1,          \
                                        colvarvalue const &x2) const    \
  {                                                                     \
    return 2.0 * (x1.real_value - x2.real_value);                       \
  }                                                                     \
                                                                        \
                                                                        \
  colvarvalue colvar::TYPE::dist2_rgrad(colvarvalue const &x1,          \
                                        colvarvalue const &x2) const    \
  {                                                                     \
    return this->dist2_lgrad(x2, x1);                                   \
  }                                                                     \


#endif
