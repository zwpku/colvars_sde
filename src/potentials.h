#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <string>
#include <vector>

class potential_function 
{
  public:
    std::string name;
    int dim;
    potential_function();
    virtual void init()=0;
    virtual void get_force(std::vector<double> &, std::vector<double> &)=0;
    virtual ~potential_function();
};

class Gaussian2d: public potential_function {
  public:
    Gaussian2d();
    virtual void init();
    virtual void get_force(std::vector<double> &, std::vector<double> &);
    virtual ~Gaussian2d();
};

class DW2d: public potential_function {
  public:
    DW2d();
    virtual void init();
    virtual void get_force(std::vector<double> &, std::vector<double> &);
    virtual ~DW2d();
};

#endif
