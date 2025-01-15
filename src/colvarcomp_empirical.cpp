// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarmodule.h"
#include "colvar.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvarcomp.h"
#include "colvarproxy.h"

colvar::empiricalcv::empiricalcv(std::string const &conf): cvc(conf) {

  empirical_cv = cvm::main()->proxy->get_empirical_cv();

  if (empirical_cv == nullptr)
  {
    cvm::error("Error: empirical cv not defined in potential: " + cvm::main()->proxy->get_pot_name() + ".\n");
  }

  set_function_type("empirical cv");

  x.type(colvarvalue::type_scalar);

  size_t n_dim = 0;
  n_dim = cvm::main()->proxy->get_dim();

  grad.resize(n_dim);
}

colvar::empiricalcv::~empiricalcv() {
}

void colvar::empiricalcv::calc_value() {
  x = empirical_cv->value(pos);
}

void colvar::empiricalcv::calc_gradients() {
}

void colvar::empiricalcv::apply_force(colvarvalue const &force) {

  empirical_cv->grad(pos, grad);

  cvm::main()->proxy->colvar_forces += grad * force.real_value ;
}

