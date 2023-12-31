#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <string>
#include <vector>

class generic_cv
{
  public:
    virtual double value(std::vector<double> &)=0;
    virtual void grad(std::vector<double> &, std::vector<double> &)=0;
};

class potential_function 
{
  public:
    potential_function();
    virtual void get_force(std::vector<double> &, std::vector<double> &)=0;
    virtual void init_state(std::vector<double> &);
    generic_cv * empirical_cv;
    ~potential_function();
    inline int dim() {
      return n_dim;
    }
  protected:
    std::string name;
    int n_dim;
};


class Gaussian2d: public potential_function {
  public:
    Gaussian2d();
    virtual void get_force(std::vector<double> &, std::vector<double> &);
    ~Gaussian2d();
    void init_state(std::vector<double> &);
};

class DW2d: public potential_function {
  public:
    DW2d();
    virtual void get_force(std::vector<double> &, std::vector<double> &);
    virtual ~DW2d();
    void init_state(std::vector<double> &);

};

class Stiff2d: public potential_function {
  public:
    Stiff2d();
    virtual void get_force(std::vector<double> &, std::vector<double> &);
    void init_state(std::vector<double> &);
    virtual ~Stiff2d();

    class cv: public generic_cv{
      public:
	double value(std::vector<double> &);
	void grad(std::vector<double> &, std::vector<double> &);
    };

  private:
    double stiff_eps;
};

class MuellerBrown: public potential_function {
  public:
    MuellerBrown();
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

#endif
