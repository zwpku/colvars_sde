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

colvar::coordinate::coordinate(std::string const &conf): cvc(conf) {
  set_function_type("coordinate");

  x.type(colvarvalue::type_scalar);

  get_keyval(conf, "index", index, 0);

  size_t n_dim = 0;
  n_dim = cvm::main()->proxy->get_dim();

  if ((index < 0) || (index >= n_dim))
  {
    cvm::error("Error: coordinate index " + cvm::to_str(index) + " is not in [0, " + cvm::to_str(n_dim-1) + "].\n");
  }

  grad.resize(n_dim);

  for (size_t i=0; i < pos.size(); i ++)
    grad[i] = 0.0;

  grad[index] = 1.0;
}

colvar::coordinate::~coordinate() {
}

void colvar::coordinate::calc_value() {

  x = pos[index];
}

void colvar::coordinate::calc_gradients() {
}

void colvar::coordinate::apply_force(colvarvalue const &force) {

  cvm::main()->proxy->colvar_forces += grad * force.real_value ;
}
