
#include "potentials.h"
#include "colvarproxy_sde.h"
#include "colvarparse.h"
#include <iostream>
#include <fstream>
#include <unistd.h>

double temp, delta_t;
int seed, n_steps, output_freq, dim;
std::string colvar_config_filename, potential_name, prefix;

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
  //parse->get_keyval(conf, "dim", dim, 2, colvarparse::parse_silent);
  parse->get_keyval(conf, "step", n_steps, 0, colvarparse::parse_silent);
  parse->get_keyval(conf, "colvars", colvar_config_filename, "", colvarparse::parse_silent);
  parse->get_keyval(conf, "potential", potential_name, "", colvarparse::parse_silent);
  parse->get_keyval(conf, "outputFreq", output_freq, -1, colvarparse::parse_silent);
  parse->get_keyval(conf, "outputPrefix", prefix, "", colvarparse::parse_silent);

  std::cout << "temp=" << temp << std::endl;
  std::cout << "delta_t=" << delta_t << std::endl;
  std::cout << "seed=" << seed << std::endl;
  std::cout << "step=" << n_steps << std::endl;
  std::cout << "colvar_config=" << colvar_config_filename << std::endl;
  std::cout << "potential=" << potential_name << std::endl;
  std::cout << "outputFreq=" << output_freq << std::endl;
  std::cout << "outputPrefix=" << prefix<< std::endl;
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
     pot_func = new Gaussian2d();
  else if (potential_name == "dw2d")
     pot_func = new DW2d();
  else if (potential_name == "mb")
     pot_func = new MuellerBrown();
  else if (potential_name == "stiff2d")
     pot_func = new Stiff2d();
  else 
  {
    std::cerr << "Error: no such potential: " << potential_name << std::endl ;
    exit(1);
  }

  dim = pot_func->dim();

  colvarproxy_sde * proxy = nullptr;

  if (!colvar_config_filename.empty())
  {
    t_inputrec *sde_inp = new t_inputrec;

    sde_inp->ref_t = temp;
    sde_inp->delta_t = delta_t;
    sde_inp->ld_seed = seed;
    sde_inp->dim = dim;
    sde_inp->empirical_cv = pot_func->empirical_cv;

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

  int i=1;

  std::ostream * output_stream = nullptr;

  if (output_freq > 0)
  {
    std::string filename;
    filename = prefix.size() > 0 ? std::string(prefix + ".txt") : std::string("traj.txt");

    // check whether already exist
    int exit_code;
    do {
      exit_code = access(filename.c_str(), F_OK);
    } while ((exit_code != 0) && (errno == EINTR));

    if (exit_code == 0) // remove the file if exists
    {
      std::rename(filename.c_str(), std::string(filename + ".BAK").c_str());
      printf("trajectory file %s removed to %s.\n", filename.c_str(), std::string(filename+".BAK").c_str());
    }

    output_stream = new std::ofstream(filename.c_str(), std::ios::binary);

    (*output_stream) << "#"+std::to_string(dim) << "\n";
  }

  pot_func->init_state(x);
  for (int step=0 ; step < n_steps; step ++)
  {
    if (step == (n_steps / 10 * i))
      {
	printf("Step=%d, %.1f%% finished.\n", step, i * 10.0);
	i++;
      }

    if ((output_freq > 0) && (step % output_freq == 0))
    {
      output_stream->setf(std::ios::scientific, std::ios::floatfield);
      (*output_stream) << step << " ";
      for (int j = 0 ; j < dim; j ++)
	(*output_stream) << x[j] << ' ';
      (*output_stream) << "\n";
    }

    pot_func->get_force(x, grad);

    for (int j=0; j < dim; j++)
    {
      r = normal_distribution(rng);
      x[j] += -1.0 * grad[j] * delta_t + coeff * r;
    }

    if (proxy)
    {
      proxy->update_data(step);
      proxy->calculateForces(x, bf);

      for (int j=0; j < dim; j++)
	x[j] += bf[j] * delta_t ;
    }
  }

  if (proxy) {
    proxy->finish();
    delete proxy;
    proxy = nullptr;
  }

  delete pot_func;
  pot_func = nullptr;

  if (output_stream)
  {
    delete output_stream;
    output_stream = nullptr;
  }

  return 0;
}
