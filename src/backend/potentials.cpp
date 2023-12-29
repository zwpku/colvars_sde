#include "potentials.h"
#include <iostream>
#include <math.h>

potential_function::potential_function() 
{
}

potential_function::~potential_function() 
{
}

void potential_function::init_state(std::vector<double> &x) 
{
  std::fill(x.begin(), x.end(), 0.0);
}

// Gaussian in x and y 
void Gaussian2d::init()
{
  name = "Gaussian 2d";
  n_dim = 2;
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
}
Gaussian2d::~Gaussian2d()
{
}

// Double-well in x, Gaussian in y 
DW2d::DW2d()
{
}

void DW2d::init()
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
}

void Stiff2d::init()
{
    name = "Stiff potential in 2d";
    n_dim = 2;
    stiff_eps = 0.5;
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
// Mueller-Brown potential

MuellerBrown::MuellerBrown()
{
}

void MuellerBrown::init()
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
