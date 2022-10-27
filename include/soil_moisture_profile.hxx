#ifndef SMP_H_INCLUDED
#define SMP_H_INCLUDED
/*
  Author: Ahmad Jan (ahmad.jan@noaa.gov)

  The code computes vertical soil moisture profiles for conceptual (e.g., CFE) and layered (e.g., LGAR) soil revervoirs.
  inputs: soil_storage, soil_storage_change_per_timestep, and soil moisture per layer when layered model is used
  output: vertical soil moisture profile for the given vertical discretization in the config file
  
  NOTE: For detailed model description please see README.md on github page
  
  @param soil_storage              [m] : soil storage (input through bmi)
  @param water_table_thickness     [m] : thickness of the water table from the bottom of the computational domain
  @param soil_moisture_profile     [-] : soil moisture content (1D vertical profile)
  @param soil_moisture_layered     [-] : soil moisture content of wetting fronts (bmi input to layered the model)
  @param soil_depths_layered       [m] : absolute depth of the wetting fronts (bmi input to layered the model)
  @param smcmax                    [-] : maximum soil moisture content (porosity)
  @param bb                        [-] : pore size distribution, beta exponent in Clapp-Hornberger (1978) function
  @param satpsi                    [m] : saturated capillary head (saturated moisture potential)
  @param ncells                    [-] : number of cells of the discretized soil column
  @param nlayers                   [-] : number of soil moisture layers, typically different than the ncells
  @param ncells_layered            [-] : number of soil wetting front, typically different than the ncells and nlayers
  @param soil_depth                [m] : depth of the computational domain
  @param last_layer_depth          [m] : depth of the last layer (for non-conceptual reservior, e.g, LGAR)
  @param soil_z                    [m] : soil discretization; 1D array of depths from the surface
  @param layers_z                  [m] : depth of layers from the surface
  @param soil_storage_model        [-] : optional models : conceptual or layered 

  @param input_var_names_model     [-] : we have different models and their inputs are different; this is to ensure that bmi inputs are consistent with the model inputs, need for ngen framework
  @param init_profile              [-] : flag for setting up initial soil moisture profile, as initially change in soil_storage is zero so we want to make sure the profile is computed at time t=0
  @param soil_storage_change_per_timestep  [m] : change in the soil storage per timestep
  @param soil_moisture_layered_option      [-] : valid for layered model only; linear or constant
  @param soil_storage_model_depth          [m] : depth of the soil storage reservoir
  @param soil_moisture_layered_bmi         [-] : if true, soil moisture content of wetting fronts is set by the bmi
 */

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <stdexcept>

using namespace std;

namespace soil_moisture_profile {

  struct soil_profile_parameters {
    int shape[3];
    double spacing[8];
    double origin[3];

    double soil_storage;
    double soil_storage_change_per_timestep;
    double water_table_thickness; // delete this
    double water_table_depth; 
    double *soil_moisture_profile;
    
    double smcmax;
    double bb;
    double satpsi;
    int ncells;
    double soil_depth;
    double last_layer_depth;
    double *soil_z;
    //double *layers_z;
   
    int soil_storage_model;
    double soil_storage_model_depth;
    int soil_moisture_layered_option;

    bool init_profile;

    // layered model
    double *soil_moisture_layered;
    double *soil_depths_layered;
    int ncells_layered;
    int max_ncells_layered;
    bool soil_depths_layered_bmi;
  };


  void SoilMoistureProfile(std::string config_file, struct soil_profile_parameters* parameters);

  void InitFromConfigFile(std::string config_file, struct soil_profile_parameters* parameters);

  // reading 1D array from the config file
  std::vector<double> ReadVectorData(std::string key);

  // update the profile for the current timestep
  void SoilMoistureProfileUpdate(struct soil_profile_parameters* parameters);

  // computes soil moisture profile for conceptual reservoir
  void SoilMoistureProfileFromConceptualReservoir(struct soil_profile_parameters* parameters);

  // computes soil moisture profile for layered-reservoir
  void SoilMoistureProfileFromLayeredReservoir(struct soil_profile_parameters* parameters);

  // computes linearly interpolated values for layered-reservoir with option = linear
  double LinearInterpolation(double z1, double z2, double t1, double t2, double z);
  
  /* 
  class SoilMoistureProfile{
  private:
    void InitializeArrays(void);
    
  public:
    int shape[3];
    double spacing[8];
    double origin[3];

    double soil_storage;
    double soil_storage_change_per_timestep;
    double water_table_thickness;
    double *soil_moisture_profile;
    double *soil_moisture_layered;

    double smcmax;
    double bb;
    double satpsi;
    int ncells;
    int nlayers;
    double soil_depth;
    double last_layer_depth;
    double *soil_z;
    double *layers_z;
   
    int soil_storage_model;
    double soil_storage_model_depth;
    std::string soil_moisture_layered_option;

    std::vector<std::string>* input_var_names_model;
    bool init_profile; 
    
    SoilMoistureProfile();
    SoilMoistureProfile(std::string config_file);

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

    // update the profile for the current timestep
    void SoilMoistureProfileUpdate();
    
    // input variable names for soil reservoir model
    std::vector<std::string>* InputVarNamesModel();
    
    ~SoilMoistureProfile();
    
  };
*/
};

#endif
