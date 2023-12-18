// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <algorithm>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarcomp.h"

#include "colvarproxy.h"


colvar::cvc::cvc()
{
  description = "uninitialized colvar component";
  b_try_scalable = true;
  sup_coeff = 1.0;
  sup_np = 1;
  period = 0.0;
  wrap_center = 0.0;
  width = 0.0;
  cvc::init_dependencies();
}


colvar::cvc::cvc(std::string const &conf)
{
  description = "uninitialized colvar component";
  b_try_scalable = true;
  sup_coeff = 1.0;
  sup_np = 1;
  period = 0.0;
  wrap_center = 0.0;
  width = 0.0;
  init_dependencies();
  colvar::cvc::init(conf);
}


int colvar::cvc::update_description()
{
  if (name.size() > 0) {
    description = "cvc " + name;
  } else {
    description = "unnamed cvc";
  }
  if (function_type.size() > 0) {
    description += " of type \"" + function_type + "\"";
  } else {
    description += " of unset type";
  }
  return COLVARS_OK;
}


int colvar::cvc::set_function_type(std::string const &type)
{
  function_type = type;
  if (function_types.size() == 0) {
    function_types.push_back(function_type);
  } else {
    if (function_types.back() != function_type) {
      function_types.push_back(function_type);
    }
  }
  update_description();

  for (size_t i = function_types.size()-1; i > 0; i--) {
    cvm::main()->cite_feature(function_types[i]+" colvar component"+
                              " (derived from "+function_types[i-1]+")");
  }
  cvm::main()->cite_feature(function_types[0]+" colvar component");
  return COLVARS_OK;
}


int colvar::cvc::init(std::string const &conf)
{
  if (cvm::debug())
    cvm::log("Initializing cvc base object.\n");

  std::string const old_name(name);

  if (name.size() > 0) {
    cvm::log("Updating configuration for component \""+name+"\"\n");
  }

  if (get_keyval(conf, "name", name, name)) {
    if ((name != old_name) && (old_name.size() > 0)) {
      cvm::error("Error: cannot rename component \""+old_name+
                 "\" after initialization (new name = \""+name+"\")",
                 COLVARS_INPUT_ERROR);
      name = old_name;
    }
  }
  update_description();

  get_keyval(conf, "componentCoeff", sup_coeff, sup_coeff);
  get_keyval(conf, "componentExp", sup_np, sup_np);
  if (sup_coeff != 1.0 || sup_np != 1) {
    cvm::main()->cite_feature("Linear and polynomial combination of colvar components");
  }
  // TODO these could be condensed into get_keyval()
  register_param("componentCoeff", reinterpret_cast<void *>(&sup_coeff));
  register_param("componentExp", reinterpret_cast<void *>(&sup_np));

  get_keyval(conf, "period", period, period);
  get_keyval(conf, "wrapAround", wrap_center, wrap_center);
  // TODO when init() is called after all constructors, check periodic flag
  register_param("period", reinterpret_cast<void *>(&period));
  register_param("wrapAround", reinterpret_cast<void *>(&wrap_center));

  get_keyval_feature(this, conf, "debugGradients",
                     f_cvc_debug_gradient, false, parse_silent);

  bool b_no_PBC = !is_enabled(f_cvc_pbc_minimum_image); // Enabled by default
  get_keyval(conf, "forceNoPBC", b_no_PBC, b_no_PBC);
  if (b_no_PBC) {
    disable(f_cvc_pbc_minimum_image);
  } else {
    enable(f_cvc_pbc_minimum_image);
  }

  // Attempt scalable calculations when in parallel? (By default yes, if available)
  get_keyval(conf, "scalable", b_try_scalable, b_try_scalable);

  if (cvm::debug())
    cvm::log("Done initializing cvc base object.\n");

  return cvm::get_error();
}


int colvar::cvc::init_total_force_params(std::string const &conf)
{
  if (cvm::get_error()) return COLVARS_ERROR;

  if (get_keyval_feature(this, conf, "oneSiteSystemForce",
                         f_cvc_one_site_total_force, is_enabled(f_cvc_one_site_total_force))) {
    cvm::log("Warning: keyword \"oneSiteSystemForce\" is deprecated: "
             "please use \"oneSiteTotalForce\" instead.\n");
  }
  if (get_keyval_feature(this, conf, "oneSiteTotalForce",
                         f_cvc_one_site_total_force, is_enabled(f_cvc_one_site_total_force))) {
    cvm::log("Computing total force on group 1 only\n");
  }

  return COLVARS_OK;
}


int colvar::cvc::init_dependencies() {
  size_t i;
  // Initialize static array once and for all
  if (features().size() == 0) {
    for (i = 0; i < colvardeps::f_cvc_ntot; i++) {
      modify_features().push_back(new feature);
    }

    init_feature(f_cvc_active, "active", f_type_dynamic);
//     The dependency below may become useful if we use dynamic atom groups

    init_feature(f_cvc_periodic, "periodic", f_type_static);

    init_feature(f_cvc_width, "defined_width", f_type_static);

    init_feature(f_cvc_lower_boundary, "defined_lower_boundary", f_type_static);

    init_feature(f_cvc_upper_boundary, "defined_upper_boundary", f_type_static);

    init_feature(f_cvc_gradient, "gradient", f_type_dynamic);

    init_feature(f_cvc_explicit_gradient, "explicit_gradient", f_type_static);

    init_feature(f_cvc_inv_gradient, "inverse_gradient", f_type_dynamic);
    require_feature_self(f_cvc_inv_gradient, f_cvc_gradient);

    init_feature(f_cvc_debug_gradient, "debug_gradient", f_type_user);
    require_feature_self(f_cvc_debug_gradient, f_cvc_gradient);
    require_feature_self(f_cvc_debug_gradient, f_cvc_explicit_gradient);

    init_feature(f_cvc_Jacobian, "Jacobian_derivative", f_type_dynamic);
    require_feature_self(f_cvc_Jacobian, f_cvc_inv_gradient);

    // Compute total force on first site only to avoid unwanted
    // coupling to other colvars (see e.g. Ciccotti et al., 2005)
    init_feature(f_cvc_one_site_total_force, "total_force_from_one_group", f_type_user);

    init_feature(f_cvc_pbc_minimum_image, "use_minimum-image_with_PBCs", f_type_user);

    // check that everything is initialized
    for (i = 0; i < colvardeps::f_cvc_ntot; i++) {
      if (is_not_set(i)) {
        cvm::error("Uninitialized feature " + cvm::to_str(i) + " in " + description);
      }
    }
  }

  // Initialize feature_states for each instance
  // default as available, not enabled
  // except dynamic features which default as unavailable
  feature_states.reserve(f_cvc_ntot);
  for (i = 0; i < colvardeps::f_cvc_ntot; i++) {
    bool avail = is_dynamic(i) ? false : true;
    feature_states.push_back(feature_state(avail, false));
  }

  // Features that are implemented by all cvcs by default
  // Each cvc specifies what other features are available
  feature_states[f_cvc_active].available = true;
  feature_states[f_cvc_gradient].available = true;

  // CVCs are enabled from the start - get disabled based on flags
  enable(f_cvc_active);

  // Explicit gradients are implemented in most CVCs. Exceptions must be specified explicitly.
  enable(f_cvc_explicit_gradient);

  // Use minimum-image distances by default
  enable(f_cvc_pbc_minimum_image);

  // Features that are implemented by default if their requirements are
  feature_states[f_cvc_one_site_total_force].available = true;

  return COLVARS_OK;
}


int colvar::cvc::setup()
{
  update_description();
  return COLVARS_OK;
}


colvar::cvc::~cvc()
{
  free_children_deps();
  remove_all_children();
}


void colvar::cvc::init_scalar_boundaries(cvm::real lb, cvm::real ub)
{
  enable(f_cvc_lower_boundary);
  lower_boundary.type(colvarvalue::type_scalar);
  lower_boundary.real_value = lb;
  enable(f_cvc_upper_boundary);
  upper_boundary.type(colvarvalue::type_scalar);
  upper_boundary.real_value = ub;
  register_param("lowerBoundary", reinterpret_cast<void *>(&lower_boundary));
  register_param("upperBoundary", reinterpret_cast<void *>(&upper_boundary));
}



colvarvalue const *colvar::cvc::get_param_grad(std::string const &param_name)
{
  colvarvalue const *ptr =
    reinterpret_cast<colvarvalue const *>(get_param_grad_ptr(param_name));
  return ptr != NULL ? ptr : NULL;
}


int colvar::cvc::set_param(std::string const &param_name,
                           void const *new_value)
{
  if (param_map.count(param_name) > 0) {

    // TODO When we can use C++11, make this a proper function map
    if (param_name.compare("componentCoeff") == 0) {
      sup_coeff = *(reinterpret_cast<cvm::real const *>(new_value));
    }
    if (param_name.compare("componentExp") == 0) {
      sup_np = *(reinterpret_cast<int const *>(new_value));
    }
    if (is_enabled(f_cvc_periodic)) {
      if (param_name.compare("period") == 0) {
        period = *(reinterpret_cast<cvm::real const *>(new_value));
      }
      if (param_name.compare("wrapAround") == 0) {
        wrap_center = *(reinterpret_cast<cvm::real const *>(new_value));
      }
    }
  }

  return colvarparams::set_param(param_name, new_value);
}


void colvar::cvc::read_data()
{
}



void colvar::cvc::calc_force_invgrads()
{
  cvm::error("Error: calculation of inverse gradients is not implemented "
             "for colvar components of type \""+function_type+"\".\n",
             COLVARS_NOT_IMPLEMENTED);
}


void colvar::cvc::calc_Jacobian_derivative()
{
  cvm::error("Error: calculation of inverse gradients is not implemented "
             "for colvar components of type \""+function_type+"\".\n",
             COLVARS_NOT_IMPLEMENTED);
}



void colvar::cvc::debug_gradients()
{
  // this function should work for any scalar cvc:
  // the only difference will be the name of the atom group (here, "group")
  // NOTE: this assumes that groups for this cvc are non-overlapping,
  // since atom coordinates are modified only within the current group

  cvm::log("Debugging gradients for " + description);

  return;
}


cvm::real colvar::cvc::dist2(colvarvalue const &x1,
                             colvarvalue const &x2) const
{
  return x1.dist2(x2);
}


colvarvalue colvar::cvc::dist2_lgrad(colvarvalue const &x1,
                                     colvarvalue const &x2) const
{
  return x1.dist2_grad(x2);
}


colvarvalue colvar::cvc::dist2_rgrad(colvarvalue const &x1,
                                     colvarvalue const &x2) const
{
  return x2.dist2_grad(x1);
}


void colvar::cvc::wrap(colvarvalue & /* x_unwrapped */) const
{
  return;
}



// Static members

std::vector<colvardeps::feature *> colvar::cvc::cvc_features;
