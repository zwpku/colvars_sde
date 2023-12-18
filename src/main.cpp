#include "colvarproxy_sde.h"
#include <iostream>

void force(std::vector<double> &x, std::vector<double> &grad)
{
  grad[0] = x[0] - 1.0;
  grad[1] = x[1];
}

int main()
{
  t_inputrec *sde_inp = new t_inputrec;

  sde_inp->ref_t = 300.0;
  sde_inp->delta_t = 0.1;
  sde_inp->ld_seed = 300;
  sde_inp->dim = 2;

  std::string prefix, filename_config("./colvars.in");

  colvarproxy_sde * proxy = new colvarproxy_sde();
  proxy->init(sde_inp, 0, prefix, filename_config);

  double delta_t = 0.1, r;
  int n_steps = 5;
  std::vector<double> x(2), grad(2), bf(2);

  // random number generation.
  std::default_random_engine rng;   
  std::normal_distribution<double> normal_distribution(0.0,1.0);

  for (int step=0 ; step < n_steps; step ++)
  {
    force(x, grad);

    r = normal_distribution(rng);
    x[0] += -1.0 * grad[0] * delta_t + r;
    x[1] += -1.0 * grad[1] * delta_t + r;

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
