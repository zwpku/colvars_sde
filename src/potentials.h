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

class Stiff2d: public potential_function {
  public:
    Stiff2d();
    virtual void init();
    virtual void get_force(std::vector<double> &, std::vector<double> &);
    virtual ~Stiff2d();
  private:
    double stiff_eps;
};

class MuellerBrown: public potential_function {
  public:
    MuellerBrown();
    virtual void init();
    virtual void get_force(std::vector<double> &, std::vector<double> &);
    virtual ~MuellerBrown();

  private:
    const double a[4]={-1, -1, -6.5, 0.7};
    const double b[4]={0, 0, 11, 0.6};
    const double c[4]={-10, -10, -6.5, 0.7};
    const double A[4]={-200, -100, -170, 15}; 
    const double xc[4]={1.0, 0, -0.5, -1.0};
    const double yc[4]={0.0, 0.5, 1.5, 1.0};
};

#endif
