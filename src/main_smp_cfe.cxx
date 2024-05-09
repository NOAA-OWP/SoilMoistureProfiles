/*
  author: Ahmad Jan Khattak
  email : ahmad.jan@noaa.gov
  date : May 2024
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <iomanip>

#include "bmi.h"
#include "cfe.h"
#include "bmi_cfe.h"

#include "../bmi/bmi.hxx"
#include "bmi_soil_moisture_profile.hxx"
#include "soil_moisture_profile.hxx"

#define SUCCESS 0


/*
  This pseudo-framework couples CFE and SoilMoistureProfiles modules to compute water table
  and soil moisture profile
*/

/*****************************************************************
 * Function to pass the parameters from CFE to SoilMoistureProfiles
 ****************************************************************/
void pass_data_from_cfe_to_smp(Bmi *cfe_bmi_model, BmiSoilMoistureProfile *smp_bmi_model) {

  double storage = 0.0;
  double storage_change = 0.0;
  double *storage_ptr = &storage;
  double *storage_change_ptr = &storage_change;
  cfe_bmi_model->get_value(cfe_bmi_model, "SOIL_STORAGE", storage_ptr);
  cfe_bmi_model->get_value(cfe_bmi_model, "SOIL_STORAGE_CHANGE", storage_change_ptr);

  smp_bmi_model->SetValue("soil_storage",storage_ptr);
  smp_bmi_model->SetValue("soil_storage_change",storage_change_ptr);

}

int main(int argc, char *argv[])
{
  BmiSoilMoistureProfile smp_bmi;

  if (argc != 3) {
    printf("Usage: ./run_smp.sh CFE \n\n");
    printf("Run soil moisture profile model through its BMI with a configuration file.\n");
    //printf("Output is written to the file `bmi_file.out`.\n");
    return SUCCESS;
  }

  /************************************************************************
   * Allocating memory to store the entire CFE BMI structure
  ************************************************************************/
  printf("\n Allocating memory to CFE BMI model structure ... \n");
  Bmi *cfe_bmi_model = (Bmi *) malloc(sizeof(Bmi));


  /************************************************************************
   * Registering the BMI model for CFE
  ************************************************************************/
  printf("Registering BMI CFE model\n");
  register_bmi_cfe(cfe_bmi_model);

  /************************************************************************
   * Initializing the BMI model for CFE and AORC and Freeze-thaw model
  ************************************************************************/
  printf("Initializeing BMI CFE %s \n", argv[1]);
  const char *cfg_file_cfe = argv[1];
  cfe_bmi_model->initialize(cfe_bmi_model, cfg_file_cfe);

  printf("Initializeing BMI SMP model %s \n", argv[2]);
  smp_bmi.Initialize(argv[2]);


  /************************************************************************
   * Get the information from the configuration here in Main
  ************************************************************************/
  printf("Get the information from the configuration here in Main\n");
  cfe_state_struct *cfe;
  cfe = (cfe_state_struct *) cfe_bmi_model->data;

  /************************************************************************
    This is the basic process for getting the SoilMoistureProfile and CFE to share data
    through BMI interface
    1. Update the CFE
    2. Get data from the CFE and Set variables in the Soil Moisture Profile
    3. Update the Soil Moisture Profile (which computes 1D soil moisture profile and water table depth)
  ************************************************************************/

  int nz = 20; // number of cells for soil discretization
  double *smc = new double[nz];
  double water_table;
  double soil_moisture_fraction;

  int nstep = 10;

  for (int i = 0; i < nstep; i++) {
      cfe_bmi_model->update(cfe_bmi_model);

      pass_data_from_cfe_to_smp(cfe_bmi_model, &smp_bmi);   // Get and Set values

      smp_bmi.Update();

      smp_bmi.GetValue("soil_moisture_profile",&smc[0]);
      smp_bmi.GetValue("soil_water_table",&water_table);
      smp_bmi.GetValue("soil_moisture_fraction",&soil_moisture_fraction);

      //for (int k = 0; k < 20; k++)
      //    std::cout <<"soil_moisture ("<< k << ") = "<< smc[k] <<"\n";
      std::cout<<"water table depth [m] = "<<water_table<<"\n";
  }

  smp_bmi.Finalize();
  cfe_bmi_model->finalize(cfe_bmi_model);

  return SUCCESS;
}
