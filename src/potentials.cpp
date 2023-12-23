#include "potentials.h"
#include <iostream>

potential_function::potential_function() 
{
}

potential_function::~potential_function() 
{
}

void Gaussian2d::init()
{
  name = "Gaussian 2d";
  dim = 2;
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


DW2d::DW2d()
{
}

void DW2d::init()
  {
    name = "Double well in x, Gaussian in y";
    dim = 2;
  }

void DW2d::get_force(std::vector<double> &x, std::vector<double> &grad)
{
  grad[0] = x[0] * (x[0]*x[0] - 1.0);
  grad[1] = x[1];
}

DW2d::~DW2d()
{
}

