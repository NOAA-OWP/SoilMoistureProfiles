#ifndef MODELCXX_C_INCLUDED
#define MODELCXX_C_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cmath> // fmin and fmax (C) or include algorithm for c++ min and max functions


#include "../include/model_cxx.hxx"


model_cxx::
ModelCXX::ModelCXX(std::string config_file)
{
  
  InitFromConfigFile(config_file);

  /*
    we can do many more things here needed at the initialization step.
    e.g., set default values, allocate/initialize arrays, call other functions
  */
  
}

/*
  Read and initialize values from configuration file
  @input  a   (double)  : coefficient of x^3
  @input  b   (double)  : coefficient of x^2
  @input  c   (double)  : coefficient of x^1
  @input  d   (double)  : coefficient of x^0
  @input  x0  (double)  : initial guess
*/


void model_cxx::
ModelCXX::InitFromConfigFile(string config_file)
{ 
  ifstream fp;
  fp.open(config_file);
  
  bool is_a_set = false;
  bool is_b_set = false;
  bool is_c_set = false;
  bool is_d_set = false;
  bool is_x0_set = false;

  
  while (fp) {

    string line;
    string param_key, param_value, param_unit;
    
    getline(fp, line);
   
    int loc_eq = line.find("=") + 1;
    int loc_u = line.find("[");
    param_key = line.substr(0,line.find("="));

    bool is_unit = line.find("[") != string::npos;

    if (is_unit)
      param_unit = line.substr(loc_u,line.find("]")+1);
    else
      param_unit = "";

    param_value = line.substr(loc_eq,loc_u - loc_eq);
    
    if (param_key == "a") {
      this->a = stod(param_value);
      is_a_set = true;
      continue;
    }
    else if (param_key == "b") {
      this->b = stod(param_value);
      is_b_set = true;
      continue;
    }
    else if (param_key == "c") {
      this->c = stod(param_value);
      is_c_set = true;
      continue;
    }
    else if (param_key == "d") {
      this->d = stod(param_value);
      is_d_set = true;
      continue;
    }
    else if (param_key == "x0") {
      this->x0 = stod(param_value);
      is_x0_set = true;
      continue;
    }
    else if (param_key == "verbosity") {
      verbosity = param_value;
      continue;
    }
  }
  
  fp.close();
  
  if (!is_a_set) {
    stringstream errMsg;
    errMsg << "coefficient \'a\' not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }
  if (!is_b_set) {
    stringstream errMsg;
    errMsg << "coefficient \'b\' not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }
  if (!is_c_set) {
    stringstream errMsg;
    errMsg << "coefficient \'c\' not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }
  if (!is_d_set) {
    stringstream errMsg;
    errMsg << "coefficient \'d\' not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }

  if (!is_x0_set) {
    printf("Initial guess not provided in the config file. Default is 0.0 \n");
  }


}


void model_cxx::
ModelCXX::Advance()
{

  FindRootCubicFunction();
  
}



void model_cxx::ModelCXX::
FindRootCubicFunction()
{
  // local variables, for convience only
  double a = this->a;
  double b = this->b;
  double c = this->c;
  double d = this->d;
  double x0 = this->x0;
  
  std::cout<<"Coefficient sets = "<< a <<" "<< b <<" "<< c <<" "<< d <<"\n";
  std::cout<<"Initial guess = "<< x0 <<"\n";
    
  double tolerance = 1.E-10;

  int count = 0;
  
  double x1 = x0; // initial guess
  double x2;
  double f, df;
  double diff;
  
  do {
    
    f  =  a * pow(x1,3.) + b * pow(x1,2.) + c * x1 + d; // the cubic function
    df = 3 * a * pow(x1,2.) + 2 * b * x1 + c;          // the derivative

    df = (df == 0.0) ? 0.0001 : df;
    
    if (df == 0.0) {
      //printf("Derivative of the function is zero. \n");
      //exit(0);
      //df = 0.00001;
    }
      
    x2 = x1 - f/df;

    diff = fabs(x2 - x1); // difference betweent the previous and new root
    x1 = x2;              // update the old root

    count++;
  } while(diff > tolerance);

  this->root = x1;
  std::cout<<"Root found = "<<this->root<<"\n"<<"Total number of iterations = "<<count<<"\n";
}

model_cxx::ModelCXX::
~ModelCXX() {}

#endif
