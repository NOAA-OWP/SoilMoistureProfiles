#ifndef MODELC_H_INCLUDED
#define MODELC_H_INCLUDED
/*
  author : Ahmad Jan (ahmad.jan@noaa.gov)
  date   : January 20, 2023
  
  Example problem:

  Find root of a cubic equation using Newton-Raphson method

  f(x) = a * x^3 + b * x^2 + c * x^1 + d
  
  @param a    : coefficient of x^3
  @param b    : coefficient of x^2
  @param c    : coefficient of x^1
  @param d    : coefficient of x^0 (constant)
  @param x0   : initial guess
  @param root : root of the equation

*/

/*
  Code : C language
  The code implements Newton-Raphson method (root finding method)
  to find a root of a cubic equation
  the namespace 'model_c' encloses a struct and functions
  - ModelC function is called from BMI (bmi_model_c) Initialize function to initialize the model's state
  - Advance function is called from BMI (bmi_model_c) to advance the model's state
  - struct 'parameters' is accessible through a pointer in BMI (bmi_model_c) to set/get BMI input/output parameters
  - remaining functions are used locally
*/

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <stdexcept>

using namespace std;

extern string verbosity;

namespace model_c {

  struct parameters {
    double a, b, c, d;
    double root;
    double x0;
  };

  // and initializes model's state at the first time step
  void ModelC(std::string config_file, struct parameters* parameters);

  void InitFromConfigFile(std::string config_file, struct parameters* parameters);

  // advance the model (from current_time to current_time + model_timestep)
  void Advance(struct parameters* parameters);

  // computes root using Netwon-Raphson method
  void FindRootCubicFunction(struct parameters* parameters);
  

};

#endif
