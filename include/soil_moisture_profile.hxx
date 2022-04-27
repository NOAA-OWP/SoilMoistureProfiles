#ifndef SMCP_H_INCLUDED
#define SMCP_H_INCLUDED
/*
  Author: Ahmad Jan (ahmad.jan@noaa.gov)

  The code computes vertical soil moisture profiles for conceptual (e.g., CFE) and layered (e.g., LGAR) soil revervoirs.
  inputs: soil_storage, soil_storage_change_per_timestep, and soil moisture per layer when layered model is used
  output: vertical soil moisture profile for the given vertical discretization in the config file
  
  NOTE: For detailed model description please see README.md on github page
  
  @param soil_storage_m // soil storage [m]
  @param soil_storage_change_per_timestep_m // change in the soil storage per timestep [m]
  @param water_table_thickness_m // thickness of the water table from the bottom of the computational domain [m]
  @param soil_moisture_profile; // soil moisture content vertical profile [-] (output) 
  @param soil_moisture_layered; // layered-soil moisture content [-], input

  @param smcmax; //maximum soil moisture content (porosity)
  @param bb;  // pore size distribution [-], beta exponent in Clapp-Hornberger (1978)
  @param satpsi; // saturated capillary head (saturated moisture potential) [m]
  @param ncells; // number of cells of the discretized soil column
  @param nlayers; // numer of soil moisture layers
  @param soil_depth; //depth of the column/domain
  @param last_layer_depth; // depth of the last layer (for non-conceptual reservior; LGAR)
  @param soilZ; // soil discretization; 1D array of depths from the surface
  @param layersZ; // depth of each layer from the surface
  @param soil_storage_model; // models : conceptual or layered 
  @param soil_moisture_profile_option; // valid for layered model only; linear or constant
  @param soil_moisture_profile_option_bmi; // option provided as an output if needed by other models that which model was used to compute the proifle.. do we need it? not sure, but it let's keep it for now 
  @param input_var_names_model; // we have different models and their inputs are different; this is to ensure that bmi inputs are consistent with the model inputs; need for ngen framework

    bool init_profile; //flag for setting up initial soil moisture profile, as initially change in soil_storage is zero so we want to make sure the profile is computed at time t=0
 */

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

using namespace std;

namespace smc_profile {
  
  class SMCProfile{
  private:
    void InitializeArrays(void);
    
  public:
    int shape[3];
    double spacing[8];
    double origin[3];

    // Model input/output variables
    double soil_storage_m; // soil storage [m]
    double soil_storage_change_per_timestep_m; // change in the soil storage per timestep [m]
    //double water_table_depth_m; // depth to water table from the surface [m]
    double water_table_thickness_m; // thickness of the water table from the bottom of the computational domain [m]
    double *soil_moisture_profile; // soil moisture content vertical profile [-] (output) 
    double *soil_moisture_layered; // layered-soil moisture content [-], input

    // model static parameters from config file
    double smcmax; //maximum soil moisture content (porosity)
    double bb;  // pore size distribution [-], beta exponent in Clapp-Hornberger (1978)
    double satpsi; // saturated capillary head (saturated moisture potential) [m]
    int ncells; // number of cells of the discretized soil column
    int nlayers; // numer of soil moisture layers
    
    double soil_depth; //depth of the column/domain
    double last_layer_depth; // depth of the last layer (for non-conceptual reservior; LGAR)
    double *soilZ; // soil discretization; 1D array of depths from the surface
    double *layersZ; // depth of each layer from the surface
   
    std::string soil_storage_model; // models : conceptual or layered 
    std::string soil_moisture_profile_option; // valid for layered model only; linear or constant
    int soil_moisture_profile_option_bmi; // option provided as an output if needed by other models that which model was used to compute the proifle.. do we need it? not sure, but it let's keep it for now 

    std::vector<std::string>* input_var_names_model; // we have different models and their inputs are different; this is to ensure that bmi inputs are consistent with the model inputs; need for ngen framework

    bool init_profile; //flag for setting up initial soil moisture profile, as initially change in soil_storage is zero so we want to make sure the profile is computed at time t=0
    
    SMCProfile();
    SMCProfile(std::string config_file);

    // initialize
    void InitFromConfigFile(std::string config_file);

    // reading 1D array from the config file
    std::vector<double> ReadVectorData(std::string key);

    // computes soil moisture profile for conceptual reservoir
    void SoilMoistureProfileFromConceptualReservoir();

    // computes soil moisture profile for layered-reservoir
    void SoilMoistureProfileFromLayeredReservoir();

    // computes linearly interpolated values for layered-reservoir with option = linear
    double LinearInterpolation(double z1, double z2, double t1, double t2, double z);

    // input variable names for soil reservoir model
    std::vector<std::string>* InputVarNamesModel();
    
    ~SMCProfile();
    
  };

};

#endif
