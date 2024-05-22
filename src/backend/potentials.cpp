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
    std::cout << "-" << config_key << ":  " << description << std::endl;
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
  x[0] = 0.0;
  x[1] = 0.0;
}

double Gaussian2d::get_potential(std::vector<double> &x)
{
  return 0.5 * x[0]*x[0] + 0.5 * x[1] * x[1];
}

void Gaussian2d::get_force(std::vector<double> &x, std::vector<double> &grad)
{
  grad[0] = x[0];
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
  a = 1.0;
  b = 3.0;
  c = 0.2;
}

void DW2d::init_state(std::vector<double> &x)
{
    x[0] = 1.0;
    x[1] = 0.0;
}

// Potential:  1/4 * (x^2-1)^2 + 1/2 * (a+b*exp(-x^2/c)) * y^2
double DW2d::get_potential(std::vector<double> &x)
{
  return 0.25 * (x[0]*x[0] - 1.0) * (x[0]*x[0] - 1.0) + 0.5 * (a + b*exp(-x[0]*x[0]/c)) * x[1]*x[1];
}

void DW2d::get_force(std::vector<double> &x, std::vector<double> &grad)
{
  grad[0] = x[0] * (x[0]*x[0] - 1.0) - b*x[1]*x[1]* x[0]/c * exp(-x[0]*x[0]/c);
  grad[1] = (a+b*exp(-x[0]*x[0]/c)) * x[1] ;
}

DW2d::~DW2d()
{
}

four_well::four_well()
{
  name = "four-well potential";
  n_dim = 2;
}

void four_well::init_state(std::vector<double> &x)
{
    x[0] = -1.0;
    x[1] = -1.0;
}

// Potential:  (x^2-1)^2 + 2 * (y^2-1)^2
double four_well::get_potential(std::vector<double> &x)
{
  return (x[0]*x[0] - 1.0) * (x[0]*x[0] - 1.0) + 2 * (x[1]*x[1] - 1.0) * (x[1]*x[1] - 1.0) ;
}

void four_well::get_force(std::vector<double> &x, std::vector<double> &grad)
{
  grad[0] = 4 * x[0] * (x[0]*x[0] - 1.0) ;
  grad[1] = 8 * x[1] * (x[1]*x[1] - 1.0) ;
}

four_well::~four_well()
{
}

// Stiff potential 2d  
//
// Potential: (x^2-1)^2 + 1/eps * (x^2+y-1)^2

Stiff2d::Stiff2d()
{
  name = "Stiff potential in 2d";
  n_dim = 2;
  stiff_eps = 0.3;
  empirical_cv = new cv();
}

void Stiff2d::init_state(std::vector<double> &x)
{
  x[0] = -1.0;
  x[1] = 0.0;
}

double Stiff2d::get_potential(std::vector<double> &x)
{
  return (x[0]*x[0]-1.0) * (x[0]*x[0]-1.0) + 1.0 / stiff_eps * (x[0]*x[0] + x[1]-1) * (x[0]*x[0]+x[1]-1);
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
  x[0] = -0.6;
  x[1] = 1.2;
}

double MuellerBrown::get_potential(std::vector<double> &x)
{
  double s = 0.0, dx, dy;

  for (int i = 0; i<4; i ++)
  {
    dx = x[0] - xc[i];
    dy = x[1] - yc[i];
    s += A[i] * exp(a[i]*dx*dx + b[i]*dx*dy + c[i]*dy*dy);
  }
  return s;
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

// Triple-well potential along circle 

TripleWellCircle::TripleWellCircle()
{
  name = "triple-well potential along circle";
  n_dim = 2;
  stiff_eps = 0.5;
}

void TripleWellCircle::init_state(std::vector<double> &x)
{
  x[0] = 1.0;
  x[1] = 0.0;
}

double TripleWellCircle::get_potential(std::vector<double> &x)
{
  double theta, r, v1;

  theta = atan2(x[1], x[0]);
  r = sqrt(x[0]*x[0]+x[1]*x[1]);

  if (theta>PI/3.0)
    v1 = pow(1-(3*theta/PI-1) * (3*theta/PI-1), 2);
  else 
    if (theta <-PI/3.0)
      v1 = pow(1-(3*theta/PI+1) * (3*theta/PI+1), 2);
    else 
      v1 = 1.0/5 * (3-2*cos(3*theta));

  return v1 + 1.0 / stiff_eps*(r-1)*(r-1)+5.0*exp(-5.0*r*r);
}

void TripleWellCircle::get_force(std::vector<double> &x, std::vector<double> &grad)
{
  double theta, r;
  double dv1_dangle, dv2_dr;

  theta = atan2(x[1], x[0]);
  r = sqrt(x[0]*x[0]+x[1]*x[1]);

   if (theta>PI/3.0)
    dv1_dangle = 12.0 / PI * (3*theta/PI-1) * ((3*theta/PI-1)*(3*theta/PI-1)- 1);
  else 
    if (theta <-PI/3.0)
      dv1_dangle = 12.0 / PI * (3*theta/PI+1) * ((3*theta/PI+1)*(3*theta/PI+1)- 1);
    else 
      dv1_dangle = 1.2 * sin(3*theta);

  dv2_dr = 2.0 * (r-1) / stiff_eps - 50 * r * exp(-5.0*r*r);

  grad[0] = -1.0 * dv1_dangle * x[1] / (r*r) + dv2_dr * x[0] / r;
  grad[1] = dv1_dangle * x[0] / (r*r) + dv2_dr * x[1] / r;
}

TripleWellCircle::~TripleWellCircle()
{
}

void define_potentials()
{
  std::cout << "\n------Available Potentials------" << std::endl;

  add_potential<Gaussian2d>("Gaussian 2d", "gaussian2d");
  add_potential<DW2d>("double-well potential in 2d", "dw2d");
  add_potential<four_well>("four-well potential in 2d", "4w");
  add_potential<Stiff2d>("stiff potential in 2d", "stiff2d");
  add_potential<MuellerBrown>("Mueller-Brown potential in 2d", "mb");
  add_potential<TripleWellCircle>("triple-well potential along circle", "triple");

  std::cout << std::endl;
}
