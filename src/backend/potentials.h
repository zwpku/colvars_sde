#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <string>
#include <vector>
#include <map>
#include <functional>

#include "colvartypes.h"

class generic_cv
{
  public:
    virtual double value(cvm::vector1d<cvm::real> &)=0;
    virtual void grad(cvm::vector1d<cvm::real> &,  cvm::vector1d<cvm::real>&)=0;
};

class potential_function 
{
  public:
    potential_function();
    virtual double get_potential(std::vector<double> &)=0;
    virtual void get_force(std::vector<double> &, std::vector<double> &)=0;
    virtual void init_state(std::vector<double> &);
    generic_cv * empirical_cv;
    ~potential_function();
    inline int dim() {
      return n_dim;
    }
    inline std::string get_name() {
      return name;
    }
  protected:
    std::string name;
    int n_dim;
};

class Gaussian2d: public potential_function {
  public:
    Gaussian2d();
    virtual double get_potential(std::vector<double> &);
    virtual void get_force(std::vector<double> &, std::vector<double> &);
    ~Gaussian2d();
    void init_state(std::vector<double> &);
};

class DW2d: public potential_function {
  public:
    DW2d();
    double a,b,c;
    virtual void get_force(std::vector<double> &, std::vector<double> &);
    virtual double get_potential(std::vector<double> &);
    virtual ~DW2d();
    void init_state(std::vector<double> &);
};

class Stiff2d: public potential_function {
  public:
    Stiff2d();
    virtual double get_potential(std::vector<double> &);
    virtual void get_force(std::vector<double> &, std::vector<double> &);
    void init_state(std::vector<double> &);
    virtual ~Stiff2d();

    class cv: public generic_cv{
      public:
	double value(cvm::vector1d<cvm::real> &);
	void grad(cvm::vector1d<cvm::real> &, cvm::vector1d<cvm::real> &);
    };

  private:
    double stiff_eps;
};

class MuellerBrown: public potential_function {
  public:
    MuellerBrown();
    virtual double get_potential(std::vector<double> &);
    virtual void get_force(std::vector<double> &, std::vector<double> &);
    void init_state(std::vector<double> &);
    virtual ~MuellerBrown();

  private:
    const double a[4]={-1, -1, -6.5, 0.7};
    const double b[4]={0, 0, 11, 0.6};
    const double c[4]={-10, -10, -6.5, 0.7};
    const double A[4]={-200, -100, -170, 15}; 
    const double xc[4]={1.0, 0, -0.5, -1.0};
    const double yc[4]={0.0, 0.5, 1.5, 1.0};
};

class TripleWellCircle: public potential_function {
  public:
    TripleWellCircle();
    virtual double get_potential(std::vector<double> &);
    virtual void get_force(std::vector<double> &, std::vector<double> &);
    void init_state(std::vector<double> &);
    virtual ~TripleWellCircle();
  private:
    double stiff_eps;
};

extern std::map<std::string, std::function<potential_function * ()>> global_potential_map;
extern std::map<std::string, std::string> global_potential_desc_map;

template <typename pot_class_name> void add_potential(char const * description, char const * config_key);
void define_potentials();

#endif
