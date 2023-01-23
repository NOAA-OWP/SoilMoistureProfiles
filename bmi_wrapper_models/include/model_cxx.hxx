#ifndef MODELCXX_H_INCLUDED
#define MODELCXX_H_INCLUDED
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
  Code : C++ language
  The code implements Newton-Raphson method (root finding method)
  to find a root of a cubic equation
  the namespace 'model_cxx' encloses the class ModelCXX
  - The constructor ModelCXX is called from BMI (bmi_model_cxx) Initialize function to initialize the model's state
  - Advance function is called from BMI (bmi_model_cxx) to advance the model's state
  - A pointer object in BMI (bmi_model_c) is used to set/get BMI input/output parameters
  - public functions are accessible in the BMI; private members are not
*/

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <stdexcept>

using namespace std;

extern string verbosity;

namespace model_cxx {

  class ModelCXX {
    
  private:
    std::string config_file;

    void InitFromConfigFile(std::string config_file);

    // computes root using Netwon-Raphson method
    void FindRootCubicFunction();
  public:
    
    double a, b, c, d;
    double root;
    double x0;

    ModelCXX() {}
    ModelCXX(std::string config_file);
    
    // advance the model (from current_time to current_time + model_timestep)
    void Advance();

    ~ModelCXX();
  };
  
};

#endif
