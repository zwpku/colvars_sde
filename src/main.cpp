
#include "colvarproxy_sde.h"
#include "colvarparse.h"
#include <iostream>
#include <fstream>

#include "potentials.h"

double temp, delta_t;
int seed, n_steps, dim;
std::string colvar_config_filename, potential_name;

void parse_input(char const  * config_filename)
{
  colvarparse * parse = new colvarparse(); // Parsing object for global options
  // open the configfile
  std::ifstream config_s(config_filename, std::ios::binary);

  // read the config file into a string
  std::string conf = "";
  std::string line;
  while (parse->read_config_line(config_s, line)) {
    // Delete lines that contain only white space after removing comments
    if (line.find_first_not_of(colvarparse::white_space) != std::string::npos)
      conf.append(line+"\n");
  }

  parse->get_keyval(conf, "temp", temp, 300.0, colvarparse::parse_silent);
  parse->get_keyval(conf, "delta_t", delta_t, 0.1, colvarparse::parse_silent);
  parse->get_keyval(conf, "seed", seed, 10, colvarparse::parse_silent);
  parse->get_keyval(conf, "dim", dim, 2, colvarparse::parse_silent);
  parse->get_keyval(conf, "step", n_steps, 0, colvarparse::parse_silent);
  parse->get_keyval(conf, "colvars", colvar_config_filename, "", colvarparse::parse_silent);
  parse->get_keyval(conf, "potential", potential_name, "", colvarparse::parse_silent);

  std::cout << "temp=" << temp << std::endl;
  std::cout << "delta_t=" << delta_t << std::endl;
  std::cout << "seed=" << seed << std::endl;
  std::cout << "step=" << n_steps << std::endl;
  std::cout << "colvar_config=" << colvar_config_filename << std::endl;
  std::cout << "potential=" << potential_name << std::endl;
}

int main(int argc, char ** argv)
{
  if (argc < 2)
  {
    std::cerr << "Error: parameter file not provided!\n" ;
    exit(1);
  }
  else 
  {
    std::cout << "Read parameters from file: " << argv[1] << "\n";
  }

  double r;
  double coeff;

  std::vector<double> x(2), grad(2), bf(2);

  parse_input(argv[1]);

  potential_function *pot_func= nullptr;

  if (potential_name == "gaussian2d")
     pot_func = new Gaussian2d;
  else if (potential_name == "dw2d")
     pot_func = new DW2d;
  else {
    std::cerr << "Error: no such potential: " << potential_name << std::endl ;
    exit(1);
  }

  colvarproxy_sde * proxy = nullptr;

  if (!colvar_config_filename.empty())
  {
    t_inputrec *sde_inp = new t_inputrec;

    sde_inp->ref_t = temp;
    sde_inp->delta_t = delta_t;
    sde_inp->ld_seed = seed;
    sde_inp->dim = dim;

    std::string prefix; 

    proxy = new colvarproxy_sde();
    proxy->init(sde_inp, 0, prefix, colvar_config_filename);
    coeff = sqrt(2 * proxy->boltzmann() * temp * delta_t);
  }
  else 
  {
    double boltzman_ = 0.001987191;
    coeff = sqrt(2 * boltzman_ * temp * delta_t);
  }

  // random number generation.
  std::default_random_engine rng;   
  std::normal_distribution<double> normal_distribution(0.0,1.0);

  x[0] = 1.0;
  x[1] = 0.0;

  int i=1;

  for (int step=0 ; step < n_steps; step ++)
  {
    if (step == (n_steps / 10 * i))
      {
	printf("Step=%d, %.1f%% finished.\n", step, i * 10.0);
	i++;
      }
    pot_func->get_force(x, grad);

    r = normal_distribution(rng);
    x[0] += -1.0 * grad[0] * delta_t + coeff * r;
    r = normal_distribution(rng);
    x[1] += -1.0 * grad[1] * delta_t + coeff * r;

    if (proxy)
    {
      proxy->update_data(step);
      proxy->calculateForces(x, bf);
      x[0] += bf[0] * delta_t ;
      x[1] += bf[1] * delta_t ;
    }
  }

  if (proxy) {
    proxy->finish();
    delete proxy;
    proxy = nullptr;
  }

  delete pot_func;
  pot_func = nullptr;

  return 0;
}
