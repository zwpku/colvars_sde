// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarcomp.h"

colvar::linearCombination::linearCombination(std::string const &conf): cvc(conf) {
    // Lookup all available sub-cvcs
    for (auto it_cv_map = colvar::get_global_cvc_map().begin(); it_cv_map != colvar::get_global_cvc_map().end(); ++it_cv_map) {
        if (key_lookup(conf, it_cv_map->first.c_str())) {
            std::vector<std::string> sub_cvc_confs;
            get_key_string_multi_value(conf, it_cv_map->first.c_str(), sub_cvc_confs);
            for (auto it_sub_cvc_conf = sub_cvc_confs.begin(); it_sub_cvc_conf != sub_cvc_confs.end(); ++it_sub_cvc_conf) {
                cv.push_back((it_cv_map->second)(*(it_sub_cvc_conf)));
            }
        }
    }
    // Sort all sub CVs by their names
    std::sort(cv.begin(), cv.end(), colvar::compare_cvc);
    for (auto it_sub_cv = cv.begin(); it_sub_cv != cv.end(); ++it_sub_cv) {
        for (auto it_atom_group = (*it_sub_cv)->atom_groups.begin(); it_atom_group != (*it_sub_cv)->atom_groups.end(); ++it_atom_group) {
            register_atom_group(*it_atom_group);
        }
    }
    // Show useful error messages and prevent crashes if no sub CVC is found
    if (cv.size() == 0) {
        cvm::error("Error: the CV " + name +
                   " expects one or more nesting components.\n");
        return;
    } else {
        x.type(cv[0]->value());
        x.reset();
    }
    use_explicit_gradients = true;
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        if (!cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            use_explicit_gradients = false;
        }
    }
    if (!use_explicit_gradients) {
        disable(f_cvc_explicit_gradient);
    }
}

cvm::real colvar::linearCombination::getPolynomialFactorOfCVGradient(size_t i_cv) const {
    cvm::real factor_polynomial = 1.0;
    if (cv[i_cv]->value().type() == colvarvalue::type_scalar) {
        factor_polynomial = cv[i_cv]->sup_coeff * cv[i_cv]->sup_np * cvm::pow(cv[i_cv]->value().real_value, cv[i_cv]->sup_np - 1);
    } else {
        factor_polynomial = cv[i_cv]->sup_coeff;
    }
    return factor_polynomial;
}

colvar::linearCombination::~linearCombination() {
    // Recall the steps we initialize the sub-CVCs:
    // 1. Lookup all sub-CVCs and then register the atom groups for sub-CVCs
    //    in their constructors;
    // 2. Iterate over all sub-CVCs, get the pointers of their atom groups
    //    groups, and register again in the parent (current) CVC.
    // That being said, the atom groups become children of the sub-CVCs at
    // first, and then become children of the parent CVC.
    // So, to destruct this class (parent CVC class), we need to remove the
    // dependencies of the atom groups to the parent CVC at first.
    remove_all_children();
    // Then we remove the dependencies of the atom groups to the sub-CVCs
    // in their destructors.
    for (auto it = cv.begin(); it != cv.end(); ++it) {
        delete (*it);
    }
    // The last step is cleaning up the list of atom groups.
    atom_groups.clear();
}

void colvar::linearCombination::calc_value() {
    x.reset();
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        cv[i_cv]->calc_value();
        colvarvalue current_cv_value(cv[i_cv]->value());
        // polynomial combination allowed
        if (current_cv_value.type() == colvarvalue::type_scalar) {
            x += cv[i_cv]->sup_coeff * (cvm::pow(current_cv_value.real_value, cv[i_cv]->sup_np));
        } else {
            x += cv[i_cv]->sup_coeff * current_cv_value;
        }
    }
}

void colvar::linearCombination::calc_gradients() {
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        cv[i_cv]->calc_gradients();
        if (cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            for (size_t j_elem = 0; j_elem < cv[i_cv]->value().size(); ++j_elem) {
                for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                    for (size_t l_atom = 0; l_atom < (cv[i_cv]->atom_groups)[k_ag]->size(); ++l_atom) {
                        (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad = factor_polynomial * (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad;
                    }
                }
            }
        }
    }
}

void colvar::linearCombination::apply_force(colvarvalue const &force) {
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        // If this CV us explicit gradients, then atomic gradients is already calculated
        // We can apply the force to atom groups directly
        if (cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                (cv[i_cv]->atom_groups)[k_ag]->apply_colvar_force(force.real_value);
            }
        } else {
            // Compute factors for polynomial combinations
            cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            colvarvalue cv_force = force.real_value * factor_polynomial;
            cv[i_cv]->apply_force(cv_force);
        }
    }
}

