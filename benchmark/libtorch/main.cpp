#include <torch/torch.h>
#include <torch/script.h>
#include <iostream>
#include <cmath>
#include <random>
#include <time.h>

using namespace std;

class Potential {

  private:
    string model_file;
    int dim;
    torch::jit::script::Module nn;

    torch::Tensor input_tensor;
    torch::Tensor nn_outputs;
    torch::Tensor input_grad;

    void load_model(string & fname);

  public:
    Potential(string & fname)
    {
      dim = 2;
      model_file = fname;
      load_model(model_file);

      auto options = torch::TensorOptions().dtype(torch::kFloat32).requires_grad(true);
      input_tensor = torch::zeros({1,(long int) dim}, options);
    }
    double value(vector<double> & x);
    void gradient(vector<double> & x, vector<double>& grad);
};

void Potential::load_model(string& fname)
{
  model_file = fname ;
  try {
    cout << "loading model from: " << fname << endl;
    nn = torch::jit::load(fname.c_str());
    nn.to(torch::kCPU);
    cout << "model loaded." << endl ;
  } catch (const std::exception & e) {
    cout << "Error: couldn't load libtorch model.\n" ;
    return ;
  }
}

double Potential::value(vector<double>& x)
{
  {
    torch::NoGradGuard no_grad;
    input_tensor[0][0] = x[0];
    input_tensor[0][1] = x[1];
  }

  std::vector<torch::jit::IValue> inputs={input_tensor};

  // evaluate the value of function
  nn_outputs = nn.forward(inputs).toTensor()[0];

  return nn_outputs.item<double>() ;
}

void Potential::gradient(vector<double>& x, vector<double>& grad)
{
  {
    torch::NoGradGuard no_grad;
    input_tensor[0][0] = x[0];
    input_tensor[0][1] = x[1];
  }

  std::vector<torch::jit::IValue> inputs={input_tensor};

  // evaluate the value of function
  nn_outputs = nn.forward(inputs).toTensor()[0];

  input_grad = torch::autograd::grad({nn_outputs}, {input_tensor})[0][0];

  for (int i=0; i < dim; i ++)
    grad[i] = input_grad[i].item<double>();
}

void generate_traj(Potential & pot, vector<double> & x0, double beta, double delta_t, int N, int seed)
{
  vector<double> x(x0), grad(2);

  default_random_engine g(seed);
  normal_distribution<double> d(0.0, 1.0);

  double coeff = sqrt(2.0 * delta_t / beta) , r;

  for (int step = 0 ; step < N; step ++)
  {
    pot.gradient(x, grad);
    r = d(g);
    x[0] += -1.0 * grad[0] * delta_t + coeff * r;
    r = d(g);
    x[1] += -1.0 * grad[1] * delta_t + coeff * r;
  }
  cout << x[0] << ',' << x[1] << ' ' << pot.value(x) << endl;
}

int main(int argc, char *argv[]) {

  string fname("model.pt");

  if (argc == 2) {
    fname = argv[1];
  }

  Potential pot(fname);

  vector<double> x0(2);
  x0[0] = -1.0; x0[1] = 0.0;

  double beta = 1.0, delta_t = 0.01;
  int N = 1000000;

  clock_t start, end;

  start = clock();
  generate_traj(pot, x0, beta, delta_t, N, 0);
  end = clock();

  double duration_sec = double (end-start) / CLOCKS_PER_SEC ; 

  cout << "Runtime: " << duration_sec << " sec." << endl ;

  return 0 ;
}

