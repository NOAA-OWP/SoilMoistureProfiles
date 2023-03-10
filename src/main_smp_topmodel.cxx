/*
  author: Ahmad Jan
  email : ahmad.jan@noaa.gov
  date : March 2023
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


#include "../topmodel/include/topmodel.h"
#include "../topmodel/include/bmi.h"
#include "../topmodel/include/bmi_topmodel.h"

#include "../bmi/bmi.hxx"
#include "../include/bmi_soil_moisture_profile.hxx"
#include "../include/soil_moisture_profile.hxx"


/*
  This pseudo-framework couples Topmodel and SoilMoistureProfiles modules to compute water table
  and soil moisture profile using two schemes based on catchment deficit, and baseflow and recharge fluxes.
  More details can be found on the main github repo, however, the two methods implemented here are based on
  Method 1 : Eq. (15) in Franchini et al. (1996))
  Method 2 : Eq. (2) in Blazkova et al. (2002)
*/

/***************************************************************
    Function to pass the parameters from Topmodel to SoilMoistureProfiles
***************************************************************/
void pass_data_from_topmodel_to_smc(Bmi *topmodel_bmi, BmiSoilMoistureProfile *smc_bmi) {
  
  
  double Qb;      // baseflow in topmodel
  double Qv;      // flow to saturated zone from unsaturated zone
  double deficit; // catchment soil moisture deficit

  std::string Qb_name = "land_surface_water__baseflow_volume_flux";
  std::string Qv_name = "soil_water_root-zone_unsat-zone_top__recharge_volume_flux";
  std::string Sbar_name = "soil_water__domain_volume_deficit";
  topmodel_bmi->get_value(topmodel_bmi, Qb_name.c_str(), &Qb);
  topmodel_bmi->get_value(topmodel_bmi, Qv_name.c_str(), &Qv);
  topmodel_bmi->get_value(topmodel_bmi, Sbar_name.c_str(), &deficit);

  smc_bmi->SetValue("Qb_topmodel",&Qb);
  smc_bmi->SetValue("Qv_topmodel",&Qv);
  smc_bmi->SetValue("global_deficit",&deficit);
  
}

/************************************************************************
    This main program is a mock framwork.
    This is not part of BMI, but acts as the driver that calls the model.
************************************************************************/
int
main(int argc, const char *argv[]) {
  
  /************************************************************************
      A configuration file is required for running this model through BMI
  ************************************************************************/
  if(argc<=1) {
    printf("make sure to include a path to config files\n");
    exit(1);
  }

  /************************************************************************
      allocating memory to store the entire BMI structure for TopModel
  ************************************************************************/
  printf("\n Allocating memory to TOPMODEL BMI model structure ... \n");
  Bmi *model = (Bmi *) malloc(sizeof(Bmi));
  
  BmiSoilMoistureProfile smp_bmi;
  
  /************************************************************************
      Registering the BMI model for Topmodel
  ************************************************************************/
  register_bmi_topmodel(model);
  printf("Registering BMI topmodel\n");
  
  /************************************************************************
      Initializing the BMI for Topmodel
  ************************************************************************/
  
  printf("Initializeing BMI Topmodel. %s \n", argv[1]);
  const char *cfg_file_topmodel = argv[1];
  model->initialize(model, cfg_file_topmodel);
  
  printf("Initializeing BMI SMP \n"); 
  const char *cfg_file_smp = argv[2];      
  smp_bmi.Initialize(cfg_file_smp); 
  
  /************************************************************************
    Get the information from the configuration here in Main
  ************************************************************************/
  printf("Get the information from the configuration here in Main\n");
  topmodel_model *topmodel;
  topmodel = (topmodel_model *) model->data;
  
  /************************************************************************
    This is the basic process for getting the SoilMoistureProfile and Topmodel to share data
    through BMI interface
    1. Update the Topmodel
    2. Get data from the Topmodel and Set variables in the Soil Moisture Profile
    3. Update the Soil Moisture Profile (which computes 1D soil moisture profile and water table depth)
  ************************************************************************/

  int nstep;
  
  if (topmodel->stand_alone == TRUE) {
    // Gather number of steps from input file
    // when in standalone mode.
    nstep = topmodel->nstep;
  }
  else {
    // Otherwise define loop here
    // Note: this is a pseudo-framework
    nstep = 720;
  }
  

  /************************************************************************
    Now loop through time and call the models with the intermediate get/set
  ************************************************************************/
  printf("looping through and calling updata\n");

  // output files -- writing water table depth, soil moisture fraction, and soil moisture profiles to separate files
  ofstream fout, fout_wt;
  fout.open("smp_data.csv");
  fout_wt.open("water_table.csv");
  fout_wt << "water_table [m]"<<",soil_moisture_fraction"<<"\n";

  int nz = 20; // number of cells for soil discretization
  double *smc = new double[nz];
  double water_table;
  double soil_moisture_fraction;
  
  for (int i = 0; i < nstep; i++) {

    model->update(model);

    pass_data_from_topmodel_to_smc(model, &smp_bmi);   // Get and Set values


    smp_bmi.Update();
    
    smp_bmi.GetValue("soil_moisture_profile",&smc[0]);
    smp_bmi.GetValue("soil_water_table",&water_table);
    smp_bmi.GetValue("soil_moisture_fraction",&soil_moisture_fraction);

    fout<<i<<",";
    for (int k = 0; k < 20; k++)
      fout <<smc[k]<<",";
    fout<< std::endl;

    fout_wt << water_table <<","<< soil_moisture_fraction <<"\n";
    //smp_bmi.PrintSoilMoistureProfile();
  }

  fout.close();
  // Run the Mass Balance check
  
  /************************************************************************
    Finalize both the Topmodel bmi
  ************************************************************************/
  printf("\n Finalizing TOPMODEL and SMP BMIs ... \n");
  model->finalize(model);

  return 0;
}

