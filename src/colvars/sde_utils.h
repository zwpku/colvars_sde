#ifndef GMX_UTILS_H
#define GMX_UTILS_H

#include "potentials.h"

struct t_inputrec 
{
    std::string pot_name;
    double delta_t;
    double ref_t;
    int64_t ld_seed;
    int dim;
    generic_cv * empirical_cv;
};

#endif
