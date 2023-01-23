#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>

#include "../bmi/bmi.hxx"

#ifdef MODELC
#include "../include/bmi_model_c.hxx"
#include "../include/model_c.hxx"
#endif


#ifdef MODELCXX
#include "../include/bmi_model_cxx.hxx"
#include "../include/model_cxx.hxx"
#endif


#define SUCCESS 0
int main(int argc, char *argv[])
{
  
#ifdef MODELC
  std::cout<<"C-based model running. \n";
  BmiModelC model;
#endif

#ifdef MODELCXX
  std::cout<<"CXX-based model running. \n";
  BmiModelCXX model;
#endif
  
  if (argc != 2) {
    printf("Usage: ./bmi_models CONFIGURATION_FILE\n\n");
    printf("BMI wrapper for C and CXX models. \n");
    printf("Output is written to the file `model_test.out`.\n");
    return SUCCESS;
  }

  FILE *fp = fopen("model_test.out", "w");
  fprintf(fp, "Configuration file = %s\n", argv[1]);
  fprintf(fp, "Initializing... ");
  
  model.Initialize(argv[1]);
  
  fprintf(fp, "done\n");
  
  {
    std::string model_name;
    model_name = model.GetComponentName();
    fprintf(fp, "%s\n", model_name.c_str());
  }
  
  {
    // inputs
    std::string var_name_a = "coefficient_a";
    std::string var_name_b = "coefficient_b";
    std::string var_name_c = "coefficient_c";
    std::string var_name_d = "coefficient_d";

    // output
    std::string var_name_root = "root";

    /****************************************************************************/
    // benchmark value(root)
    double root_check = 1.;

    // Example 1: Set values for equation : 3 * x^3 - 16 * x^2 + 23 * x - 6 = 0 (roots=3,2,1/3)
    
    double a = 3;
    double b = -16;
    double c = 23;
    double d = -6;
    
    // Example 2: Set values for equation : 2 * x^3 - 11 * x^2 + 17 * x - 6 = 0 (roots=2,3,0.5)
    /*
    double a = 2;
    double b = -11;
    double c = 17;
    double d = -6;
    */
    
    double *a_ptr = &a;
    double *b_ptr = &b;
    double *c_ptr = &c;
    double *d_ptr = &d;

    model.SetValue(var_name_a,a_ptr);
    model.SetValue(var_name_b,b_ptr);
    model.SetValue(var_name_c,c_ptr);
    model.SetValue(var_name_d,d_ptr);

    
    double root_computed;
    model.GetValue(var_name_root, &root_computed);
    
    model.Update();

    model.GetValue(var_name_root, &root_computed);
    
    std::cout<<"root: "<<root_computed<<"\n";
  }
  
  fprintf(fp, "Finalizing... ");

  model.Finalize();
  
  fprintf(fp, "done\n");
  fclose(fp);
  return SUCCESS;
}
