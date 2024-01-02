#include "potentials.h"
#include <iostream>
#include <math.h>

std::map<std::string, std::function<potential_function * ()>> global_potential_map = std::map<std::string, std::function<potential_function * ()>>() ;
std::map<std::string, std::string> global_potential_desc_map = std::map<std::string, std::string>() ;

template <typename pot_class_name>
void add_potential(char const * description, char const * config_key)
{
  if (global_potential_map.count(config_key) == 0) {
    global_potential_map[config_key] = []() {
      return new pot_class_name();
    };
    global_potential_desc_map[config_key] = std::string(description);
    std::cout << config_key << " " << description << std::endl;
  }
}

potential_function::potential_function() 
{
  empirical_cv = nullptr;
}

potential_function::~potential_function() 
{
  if (empirical_cv)
  {
    delete empirical_cv;
    empirical_cv = nullptr;
  }
}

void potential_function::init_state(std::vector<double> &x) 
{
  std::fill(x.begin(), x.end(), 0.0);
}

void Gaussian2d::init_state(std::vector<double> &x)
{
  x[0] = 1.0;
  x[1] = 0.0;
}

void Gaussian2d::get_force (std::vector<double> &x, std::vector<double> &grad)
{
  grad[0] = x[0] - 1.0;
  grad[1] = x[1];
}


Gaussian2d::Gaussian2d()
{
  name = "Gaussian 2d";
  n_dim = 2;
}

Gaussian2d::~Gaussian2d()
{
}

// Double-well in x, Gaussian in y 
DW2d::DW2d()
{
  name = "Double well in x, Gaussian in y";
  n_dim = 2;
}

void DW2d::init_state(std::vector<double> &x)
{
    x[0] = 1.0;
    x[0] = 0.0;
}

void DW2d::get_force(std::vector<double> &x, std::vector<double> &grad)
{
  grad[0] = x[0] * (x[0]*x[0] - 1.0);
  grad[1] = x[1];
}

DW2d::~DW2d()
{
}

// Stiff potential 2d  

Stiff2d::Stiff2d()
{
  name = "Stiff potential in 2d";
  n_dim = 2;
  stiff_eps = 0.5;
  empirical_cv = new cv();
}

void Stiff2d::init_state(std::vector<double> &x)
{
  x[0] = -1.0;
  x[1] = 0.0;
}

void Stiff2d::get_force(std::vector<double> &x, std::vector<double> &grad)
{
  grad[0] = 4 * x[0] * (x[0]*x[0] - 1.0 + 1.0 / stiff_eps * (x[0]*x[0] + x[1] -1));
  grad[1] = 2.0 / stiff_eps * (x[0]*x[0] + x[1] -1);
}

Stiff2d::~Stiff2d()
{
}

double Stiff2d::cv::value(cvm::vector1d<cvm::real> &x)
{
  return 0;
}

void Stiff2d::cv::grad(cvm::vector1d<cvm::real> &x, cvm::vector1d<cvm::real> & grad)
{
}

// Mueller-Brown potential

MuellerBrown::MuellerBrown()
{
  name = "Mueller-Brown potential in 2d";
  n_dim = 2;
}

void MuellerBrown::init_state(std::vector<double> &x)
{
  x[0] = -1.0;
  x[1] = 1.0;
}

void MuellerBrown::get_force(std::vector<double> &x, std::vector<double> &grad)
{
  double dx, dy;

  grad[0] = grad[1] = 0;
  for (int i = 0; i<4; i ++)
  {
    dx = x[0] - xc[i];
    dy = x[1] - yc[i];
    grad[0] += A[i] * (2.0 * a[i] * dx + b[i] * dy) * exp(a[i]*dx*dx + b[i]*dx*dy + c[i]*dy*dy);
    grad[1] += A[i] * (b[i] * dx + 2.0 * c[i] * dy) * exp(a[i]*dx*dx + b[i]*dx*dy + c[i]*dy*dy);
  }
}

MuellerBrown::~MuellerBrown()
{
}

void define_potentials()
{
  add_potential<Gaussian2d>("Gaussian 2d", "gaussian2d");
  add_potential<DW2d>("double-well potential in 2d", "dw2d");
  add_potential<Stiff2d>("stiff potential in 2d", "stiff2d");
  add_potential<MuellerBrown>("Mueller-Brown potential in 2d", "mb");
}
