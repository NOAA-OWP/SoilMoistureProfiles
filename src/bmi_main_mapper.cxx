#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "../../bmi/bmi.hxx"
#include "../include/bmi_mapper.hxx"
#include "../include/smc_mapping.hxx"

#define SUCCESS 0
int main(int argc, char *argv[])
{
  BmiMapper model;
  
  if (argc != 2) {
    printf("Usage: run_bmifrozensoilcxx CONFIGURATION_FILE\n\n");
    printf("Run the frozensoilcxx model through its BMI with a configuration file.\n");
    printf("Output is written to the file `bmifrozensoilcxx.out`.\n");
    return SUCCESS;
  }

  FILE *fp = fopen("bmi_file.out", "w");
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
    std::string var_name_id = "grid_gid_unique";
    std::string var_name_smc = "grid_SMC";
    
    //    int grid, rank, *shape;
    

    fprintf(fp, "variable = %s\n", var_name_id.c_str());
    fprintf(fp, "variable = %s\n", var_name_smc.c_str());
    
    //grid = model.GetVarGrid(var_name_s);

    //rank = model.GetGridRank(grid);
    //fprintf(fp, "rank = %d\n", rank);
    //shape = new int[rank];
    //model.GetGridShape(grid, shape);

    //fprintf(fp, "shape = %d x %d x %d\n", shape[0],1,1);

    // Set values
    int ngrids = 575;
    int *var_id = new int[ngrids];
    double *var_smc = new double[ngrids];
    
    //    model.SetValue(var_name_id,&var_id[0]);
    
    var_id = (int *)model.GetValuePtr(var_name_id);
    var_smc = (double *)model.GetValuePtr(var_name_smc);

    std::vector<double> smc_vec(ngrids);
    for (int i=0; i < ngrids; i++)
      smc_vec[i] = var_smc[i];

    double smc_min_sim = *std::min_element(smc_vec.begin(), smc_vec.end());
    double smc_max_sim = *std::max_element(smc_vec.begin(), smc_vec.end());
    double smc_mean_sim = std::accumulate(smc_vec.begin(), smc_vec.end(),0.0);
    smc_mean_sim /=smc_vec.size();

    
    //model.Update();
    // values from python for comparison
    double smc_min = 0.1040445;
    double smc_max = 0.30984;
    double smc_mean = 0.233984;
    std::cout<<"Max: "<<smc_max<<" "<<smc_max_sim<<"\n";
    std::cout<<"Min: "<<smc_min<<" "<<smc_min_sim<<"\n";
    std::cout<<"Mean: "<<smc_mean<<" "<<smc_mean_sim<<"\n";
    for (int i=0; i < 575; i++) {
      //      assert (abs(var_smc[i] - SMCT[i]) < 1.E-6);     
      fprintf(fp, "%d,%lf", var_id[i], var_smc[i]);
      fprintf(fp, "\n");
    }
    
  }
  
  fprintf(fp, "Finalizing... ");
 
  model.Finalize();
  fprintf(fp, "done\n");
  fclose(fp);
  return SUCCESS;
  
}
