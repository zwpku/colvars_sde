#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <time.h>

using namespace std;

class StiffPot {
  double eps ;

  public:
    StiffPot()
    {
      eps = 0.5;
    }
    double V(vector<double> & x)
    {
      double tmp ;
      tmp = x[0]*x[0] - 1.0;
      return tmp * tmp + 1.0 / eps * (tmp + x[1]) * (tmp + x[1]);
    }
    void gradV (vector<double> & x, vector<double>& grad)
    {
      double tmp ;
      tmp = x[0]*x[0] - 1.0;
      grad[0] = 4.0 * x[0] * (tmp + 1.0 / eps * (tmp + x[1]));
      grad[1] = 2.0 / eps * (tmp + x[1]);
    }
};

void generate_traj(vector<double> & x0, double beta, double delta_t, int N, int seed)
{
  vector<double> x(x0), grad(2);

  default_random_engine g(seed);
  normal_distribution<double> d(0.0, 1.0);

  StiffPot pot;

  double coeff = sqrt(2.0 * delta_t / beta) , r;

  for (int step = 0 ; step < N; step ++)
  {
    pot.gradV(x, grad);
    r = d(g);
    x[0] += -1.0 * grad[0] * delta_t + coeff * r;
    r = d(g);
    x[1] += -1.0 * grad[1] * delta_t + coeff * r;
  }
  cout << x[0] << ',' << x[1] << ' ' << pot.V(x) << endl;
}

int main()
{
  vector<double> x0(2);
  x0[0] = -1.0; x0[1] = 0.0;

  double beta = 1.0, delta_t = 0.01;
  int N = 1000000;

  clock_t start, end;

  start = clock();
  generate_traj(x0, beta, delta_t, N, 0);
  end = clock();

  double duration_sec = double (end-start) / CLOCKS_PER_SEC ; 

  cout << "Runtime: " << duration_sec << " sec." << endl ;

  return 0;
}
