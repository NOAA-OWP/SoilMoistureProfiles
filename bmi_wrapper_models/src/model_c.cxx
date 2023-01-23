#ifndef MODELC_C_INCLUDED
#define MODELC_C_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cmath> // fmin and fmax (C)


#include "../include/model_c.hxx"


void model_c::
ModelC(string config_file, struct parameters* parameters)
{
  
  InitFromConfigFile(config_file, parameters);

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
void model_c::
InitFromConfigFile(string config_file, struct parameters* parameters)
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
      parameters->a = stod(param_value);
      is_a_set = true;
      continue;
    }
    else if (param_key == "b") {
      parameters->b = stod(param_value);
      is_b_set = true;
      continue;
    }
    else if (param_key == "c") {
      parameters->c = stod(param_value);
      is_c_set = true;
      continue;
    }
    else if (param_key == "d") {
      parameters->d = stod(param_value);
      is_d_set = true;
      continue;
    }
    else if (param_key == "x0") {
      parameters->x0 = stod(param_value);
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


void model_c::
Advance(struct parameters* parameters)
{

  FindRootCubicFunction(parameters);
  
}



void model_c::
FindRootCubicFunction(struct parameters* parameters)
{
  // local variables, for convience only
  double a = parameters->a;
  double b = parameters->b;
  double c = parameters->c;
  double d = parameters->d;
  double x0 = parameters->x0;
  
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

  parameters->root = x1;
  
  std::cout<<"Root found = "<< parameters->root <<"\n"<<"Total number of iterations = "<< count <<"\n";

}


#endif
