#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <iomanip>

#include "../bmi/bmi.hxx"
#include "../include/bmi_soil_moisture_profile.hxx"
#include "../include/soil_moisture_profile.hxx"

#define SUCCESS 0
int main(int argc, char *argv[])
{
  BmiSoilMoistureProfile model;
  
  if (argc != 2) {
    printf("Usage: ./make_run_standalone.sh \n\n");
    printf("Run soil moisture profile model through its BMI with a configuration file.\n");
    printf("Output is written to the file `bmi_file.out`.\n");
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
    std::string var_name_s = "soil_storage";
    std::string var_name_sc = "soil_storage_change";
    std::string var_name_wt = "soil_water_table";
    std::string var_name_smc = "soil_moisture_profile";
    std::string var_name_smcl = "soil_moisture_wetting_fronts";
    std::string var_name_smc_bmi = "soil_storage_model";
    
    int grid, rank, *shape;
    double *var_s = NULL;
    double *var_sc = NULL;
    double *water_table_thickness_bmi = NULL;

    fprintf(fp, "variable = %s\n", var_name_s.c_str());
    fprintf(fp, "variable = %s\n", var_name_sc.c_str());
    fprintf(fp, "variable = %s\n", var_name_wt.c_str());
    fprintf(fp, "variable = %s\n", var_name_smc.c_str());
    fprintf(fp, "variable = %s\n", var_name_smcl.c_str());
    fprintf(fp, "variable = %s\n", var_name_smc_bmi.c_str());
    
    grid = model.GetVarGrid(var_name_smc);

    rank = model.GetGridRank(grid);
    fprintf(fp, "rank = %d\n", rank);
    shape = new int[rank];
    model.GetGridShape(grid, shape);

    fprintf(fp, "shape = %d x %d x %d\n", shape[0],1,1);

    /****************************************************************************/
    // unit test data for conceptual soil reservoir
    double soil_moisture_profile[] = {0.3375956134,0.3423608214,0.3475802559,0.3533405627,0.3597547738,0.3669739833,
				      0.3752061855,0.3847482243,0.3960434405,0.4097942568,0.4272064334,0.4390000000,
				      0.4390000000,0.4390000000,0.4390000000,0.4390000000,0.4390000000,0.4390000000,
				      0.4390000000,0.4390000000};
    double water_table_thickness = 1.50956; // in meters
    enum option { Conceptual = 1, Layered = 2};
    /****************************************************************************/
    // Set values
    double storage_m = 0.8;
    double storage_change_m = -0.000472;
    double *storage_m_ptr = &storage_m;
    double *storage_change_m_ptr = &storage_change_m;
    double smc_layers[] = {0.25, 0.15, 0.1, 0.12};

    int soil_moisture_profile_option;

    model.GetValue(var_name_smc_bmi,&soil_moisture_profile_option);

    model.SetValue(var_name_s,storage_m_ptr);

    model.SetValue(var_name_sc,storage_change_m_ptr);

    model.SetValue(var_name_smcl,&smc_layers[0]);
    
    var_s = (double *)model.GetValuePtr(var_name_s);
    var_sc = (double *)model.GetValuePtr(var_name_sc);

    std::cout<<"soil_storage_model: "<<soil_moisture_profile_option<<"\n";
    std::cout<<"storage: "<<*var_s<<"\n";
    std::cout<<"storage change: "<<*var_sc<<"\n";
    
    model.Update();

    water_table_thickness_bmi = (double *)model.GetValuePtr(var_name_wt);
    
    if (soil_moisture_profile_option == Conceptual) {
      std::cout<<"Check: water table thickness = "<<water_table_thickness<<" | water table thickness bmi = "<<*water_table_thickness_bmi<<"\n";
      double diff = std::fabs(water_table_thickness - *water_table_thickness_bmi);
      assert (diff < 1e-4);
    }
    
    // Get values
    double *var_smc = new double[shape[0]];
    
    model.GetValue(var_name_smc,&var_smc[0]);
      
    if (soil_moisture_profile_option == Conceptual) {
      std::cout<<"Referance value | Simulated value | Difference \n";
      for (int i=0; i < shape[0]; i++) {
	std::cout<< left << setw(18) << var_smc[i]
		 << setw(18) << soil_moisture_profile[i]
		 << setw(1) << abs(var_smc[i] - soil_moisture_profile[i])<<"\n";
	assert (abs(var_smc[i] - soil_moisture_profile[i]) < 1.E-6);     
	fprintf(fp, "%6.4e", var_smc[i]);
	fprintf(fp, "\n");
      }
    }
  }
  
  fprintf(fp, "Finalizing... ");

  model.Finalize();
  fprintf(fp, "done\n");
  fclose(fp);
  return SUCCESS;
}
