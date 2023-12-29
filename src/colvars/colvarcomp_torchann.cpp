// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifdef TORCH

#include "colvarmodule.h"
#include "colvar.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvarcomp.h"
#include "colvarproxy.h"

colvar::torchANN::torchANN(std::string const &conf): cvc(conf) {
  set_function_type("torchANN");

  x.type(colvarvalue::type_scalar);

  if (period != 0.0) {
    enable(f_cvc_periodic);
  }

  if ((wrap_center != 0.0) && !is_enabled(f_cvc_periodic)) {
    cvm::error("Error: wrapAround was defined in a torchANN component,"
                " but its period has not been set.\n");
    return;
  }

  std::string model_file ;
  get_keyval(conf, "modelFile", model_file, std::string(""));
  try {
    nn = torch::jit::load(model_file);
    nn.to(torch::kCPU);
    cvm::log("torch model loaded.") ;
  } catch (const std::exception & e) {
    cvm::error("Error: couldn't load libtorch model (see below).\n" + cvm::to_str(e.what()));
    return;
  }
  get_keyval(conf, "m_output_index", m_output_index, 0);
  get_keyval(conf, "doubleInputTensor", use_double_input, false);
  //get_keyval(conf, "useGPU", use_gpu, false);

  //cvm::log("Input dimension of model: " + cvm::to_str(num_inputs));

  // initialize the input tensor 
  auto options = torch::TensorOptions().dtype(torch::kFloat32).requires_grad(true);
  
  if (use_double_input) {  // set type to double
    options = options.dtype(torch::kFloat64);
    nn.to(torch::kFloat64);
    cvm::log("Model's dtype: kFloat64.");
  } else {
    cvm::log("Model's dtype: kFloat32.");
  }

  int num_inputs = 0;
  num_inputs = cvm::main()->proxy->get_dim();
  grad.resize(num_inputs);
  input_tensor = torch::zeros({1,(long int) num_inputs}, options);

  try { // test the model 
    std::vector<torch::jit::IValue> inputs={input_tensor};
    nn_outputs = nn.forward(inputs).toTensor()[0][m_output_index];
    cvm::log("Evaluating model with zero tensor succeeded.");
  } catch (const std::exception & e) {
    cvm::error("Error: evaluating model with zero tensor failed (see below).\n" + cvm::to_str(e.what()));
    return;
  } 
}

colvar::torchANN::~torchANN() {
}

void colvar::torchANN::calc_value() {

  {
    torch::NoGradGuard no_grad;
    for (size_t i= 0; i < pos.size(); ++i) 
	input_tensor[0][i] = pos[i];
  }

  std::vector<torch::jit::IValue> inputs={input_tensor};

  // evaluate the value of function
  nn_outputs = nn.forward(inputs).toTensor()[0][m_output_index];

  torch::Tensor input_grad = torch::autograd::grad({nn_outputs}, {input_tensor})[0][0];

  for (size_t i=0; i < pos.size(); i ++)
    grad[i] = input_grad[i].item<double>();

  x = nn_outputs.item<double>() ;

  this->wrap(x);

}

void colvar::torchANN::calc_gradients() {
}

void colvar::torchANN::apply_force(colvarvalue const &force) {

  cvm::main()->proxy->colvar_forces += grad * force.real_value ;
}

// Nearest-image convection is handled in the same way as in colvar::distance_z
cvm::real colvar::torchANN::dist2(colvarvalue const &x1, colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  if (is_enabled(f_cvc_periodic)) {
    cvm::real shift = cvm::floor(diff/period + 0.5);
    diff -= shift * period;
  }
  return diff * diff;
}

colvarvalue colvar::torchANN::dist2_lgrad(colvarvalue const &x1,
                                          colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  if (is_enabled(f_cvc_periodic)) {
    cvm::real shift = cvm::floor(diff/period + 0.5);
    diff -= shift * period;
  }
  return 2.0 * diff;
}

colvarvalue colvar::torchANN::dist2_rgrad(colvarvalue const &x1,
                                          colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  if (is_enabled(f_cvc_periodic)) {
    cvm::real shift = cvm::floor(diff/period + 0.5);
    diff -= shift * period;
  }
  return (-2.0) * diff;
}


void colvar::torchANN::wrap(colvarvalue &x_unwrapped) const
{
  if (!is_enabled(f_cvc_periodic)) {
    return;
  }
  cvm::real shift =
    cvm::floor((x_unwrapped.real_value - wrap_center) / period + 0.5);
  x_unwrapped.real_value -= shift * period;
}

#endif
