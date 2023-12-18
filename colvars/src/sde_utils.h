#ifndef GMX_UTILS_H
#define GMX_UTILS_H

struct t_inputrec 
{
    double delta_t;
    double ref_t;
    int64_t ld_seed;
    int dim;
};

class Force 
{
  double f;
};

class ForceProviderInput
{
public:
    /*! \brief Constructor assembles all necessary force provider input data
     *
     * \param[in]  x        Atomic positions
     * \param[in]  time     The current time in the simulation
     */
    ForceProviderInput(double x,
                       double               time
                       ) :
        x_(x),
        t_(time)
    {}

    double  x_;       //!< The atomic positions
    double               t_;       //!< The current time in the simulation
};

class ForceProviderOutput
{
public:
    /*! \brief Constructor assembles all necessary force provider output data
     *
     * \param[in,out]  forceWithVirial  Container for force and virial
     * \param[in,out]  enerd            Structure containing energy data
     */
    ForceProviderOutput(Force* force, double enerd) :
        force_(*force),
        enerd_(enerd)
    {
    }

    Force& force_; //!< Container for force 
    double enerd_;           //!< Structure containing energy data
};

#endif
