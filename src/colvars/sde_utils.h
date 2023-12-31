#ifndef GMX_UTILS_H
#define GMX_UTILS_H

#include "potentials.h"

struct t_inputrec 
{
    double delta_t;
    double ref_t;
    int64_t ld_seed;
    int dim;
    generic_cv * empirical_cv;
};

#endif
