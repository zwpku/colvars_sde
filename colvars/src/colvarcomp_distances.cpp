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
#include "colvar_rotation_derivative.h"


colvar::distance::distance(std::string const &conf)
  : cvc(conf)
{
  set_function_type("distance");
  init_as_distance();

  provide(f_cvc_inv_gradient);
  provide(f_cvc_Jacobian);
  enable(f_cvc_com_based);

  group1 = parse_group(conf, "group1");
  group2 = parse_group(conf, "group2");

  init_total_force_params(conf);
}


colvar::distance::distance()
  : cvc()
{
  set_function_type("distance");
  init_as_distance();

  provide(f_cvc_inv_gradient);
  provide(f_cvc_Jacobian);
  enable(f_cvc_com_based);
}


void colvar::distance::calc_value()
{
  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    dist_v = group2->center_of_mass() - group1->center_of_mass();
  } else {
    dist_v = cvm::position_distance(group1->center_of_mass(),
                                    group2->center_of_mass());
  }
  x.real_value = dist_v.norm();
}


void colvar::distance::calc_gradients()
{
  cvm::rvector const u = dist_v.unit();
  group1->set_weighted_gradient(-1.0 * u);
  group2->set_weighted_gradient(       u);
}


void colvar::distance::calc_force_invgrads()
{
  group1->read_total_forces();
  if (is_enabled(f_cvc_one_site_total_force)) {
    ft.real_value = -1.0 * (group1->total_force() * dist_v.unit());
  } else {
    group2->read_total_forces();
    ft.real_value = 0.5 * ((group2->total_force() - group1->total_force()) * dist_v.unit());
  }
}


void colvar::distance::calc_Jacobian_derivative()
{
  jd.real_value = x.real_value ? (2.0 / x.real_value) : 0.0;
}


void colvar::distance::apply_force(colvarvalue const &force)
{
  if (!group1->noforce)
    group1->apply_colvar_force(force.real_value);

  if (!group2->noforce)
    group2->apply_colvar_force(force.real_value);
}


simple_scalar_dist_functions(distance)


colvar::cartesian::cartesian(std::string const &conf)
  : cvc(conf)
{
  set_function_type("cartesian");

  atoms = parse_group(conf, "atoms");

  bool use_x, use_y, use_z;
  get_keyval(conf, "useX", use_x, true);
  get_keyval(conf, "useY", use_y, true);
  get_keyval(conf, "useZ", use_z, true);

  axes.clear();
  if (use_x) axes.push_back(0);
  if (use_y) axes.push_back(1);
  if (use_z) axes.push_back(2);

  if (axes.size() == 0) {
    cvm::error("Error: a \"cartesian\" component was defined with all three axes disabled.\n");
    return;
  }

  x.type(colvarvalue::type_vector);
  disable(f_cvc_explicit_gradient);
  // Don't try to access atoms if creation of the atom group failed
  if (atoms != NULL) x.vector1d_value.resize(atoms->size() * axes.size());
}


void colvar::cartesian::calc_value()
{
  size_t const dim = axes.size();
  size_t ia, j;
  for (ia = 0; ia < atoms->size(); ia++) {
    for (j = 0; j < dim; j++) {
      x.vector1d_value[dim*ia + j] = (*atoms)[ia].pos[axes[j]];
    }
  }
}


void colvar::cartesian::calc_gradients()
{
  // we're not using the "grad" member of each
  // atom object, because it only can represent the gradient of a
  // scalar colvar
}


void colvar::cartesian::apply_force(colvarvalue const &force)
{
  size_t const dim = axes.size();
  size_t ia, j;
  if (!atoms->noforce) {
    cvm::rvector f;
    for (ia = 0; ia < atoms->size(); ia++) {
      for (j = 0; j < dim; j++) {
        f[axes[j]] = force.vector1d_value[dim*ia + j];
      }
      (*atoms)[ia].apply_force(f);
    }
  }
}
