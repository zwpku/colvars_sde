#include <torch/torch.h>
#include <torch/script.h>
#include <iostream>

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
    nn = torch::jit::load(fname.c_str());
    nn.to(torch::kCPU);
    cout << "torch model loaded." << endl ;
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


int main() {

  string fname("model.pt");
  Potential pot(fname);

  vector<double> x(2), grad(2);
  x[0] = -1.0;
  x[1] = 1.0;

  cout << "potential=" << pot.value(x) << endl;

  pot.gradient(x, grad);

  cout << "gradient=" << grad[0] << ',' << grad[1] << endl;

  return 0 ;
}


