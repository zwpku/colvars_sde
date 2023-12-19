#include "colvarproxy_sde.h"
#include <iostream>

void force(std::vector<double> &x, std::vector<double> &grad)
{
  grad[0] = x[0] - 1.0;
  grad[1] = x[1];
}

int main()
{
  double delta_t = 0.01, r;
  int n_steps = 20000;
  std::vector<double> x(2), grad(2), bf(2);

  t_inputrec *sde_inp = new t_inputrec;

  sde_inp->ref_t = 300.0;
  sde_inp->delta_t = delta_t;
  sde_inp->ld_seed = 300;
  sde_inp->dim = 2;

  std::string prefix, filename_config("./colvars.in");

  colvarproxy_sde * proxy = new colvarproxy_sde();
  proxy->init(sde_inp, 0, prefix, filename_config);

  double coeff = sqrt(2 * proxy->boltzmann() * sde_inp->ref_t * delta_t);

  // random number generation.
  std::default_random_engine rng;   
  std::normal_distribution<double> normal_distribution(0.0,1.0);

  for (int step=0 ; step < n_steps; step ++)
  {
    force(x, grad);

    r = normal_distribution(rng);
    x[0] += -1.0 * grad[0] * delta_t + coeff * r;
    r = normal_distribution(rng);
    x[1] += -1.0 * grad[1] * delta_t + coeff * r;

    proxy->update_data(step);
    proxy->calculateForces(x, bf);

    x[0] += bf[0] * delta_t ;
    x[1] += bf[1] * delta_t ;
  }

  proxy->finish();
  delete proxy;
  proxy = nullptr;

  return 0;
}
