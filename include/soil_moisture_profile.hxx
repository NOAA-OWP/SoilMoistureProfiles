#ifndef SMP_H_INCLUDED
#define SMP_H_INCLUDED
/*
  Author: Ahmad Jan (ahmad.jan@noaa.gov)

  The code computes vertical soil moisture profiles for conceptual (e.g., CFE) and layered (e.g., LGAR) soil revervoirs.
  inputs: soil_storage, soil_storage_change_per_timestep, and soil moisture per layer when layered model is used
  output: vertical soil moisture profile for the given vertical discretization in the config file

  The code has also been extended to compute watertable and then soil moisture profile from soil deficit and
  fluxes (baseflow and recharge rate) in the topmodel. We call them deficit-based and iteration-based methods
  
  NOTE: For detailed model description please see README.md on github page
  
  @param soil_storage                  [m] : soil storage (input through bmi)
  @param water_table_depth             [m] : depth from the surface to the water table location
  @param soil_moisture_profile         [-] : soil moisture content (1D vertical profile)
  @param soil_moisture_wetting_fronts  [-] : soil moisture content of wetting fronts (bmi input to layered the model)
  @param soil_depth_wetting_fronts     [m] : absolute depth of the wetting fronts (bmi input to layered the model)
  @param smcmax                    [-] : maximum soil moisture content (porosity)
  @param b                         [-] : pore size distribution, beta exponent in Clapp-Hornberger (1978) function
  @param satpsi                    [m] : saturated capillary head (saturated moisture potential)
  @param ncells                    [-] : number of cells of the discretized soil column
  @param nlayers                   [-] : number of soil moisture layers, typically different than the ncells
  @param ncells_layered            [-] : number of soil wetting front, typically different than the ncells and nlayers
  @param soil_depth                [m] : depth of the computational domain
  @param soil_depth_NWM            [m] : National Water Model 3.0 soil column depth (set to 2.0 m)
  @param last_layer_depth          [m] : depth of the last layer (for non-conceptual reservior, e.g, LGAR)
  @param soil_z                    [m] : soil discretization; 1D array of depths from the surface
  @param layers_z                  [m] : depth of layers from the surface
  @param soil_storage_model        [-] : optional models : conceptual or layered

  @param input_var_names_model     [-] : we have different models and their inputs are different; this is to ensure
                                         that bmi inputs are consistent with the model inputs, need for ngen framework
  @param init_profile              [-] : flag for setting up initial soil moisture profile, as initially change in
                                         soil_storage is zero so we want to make sure the profile is computed at time t=0
  
  @param soil_storage_change_per_timestep  [m] : change in the soil storage per timestep
  @param soil_moisture_layered_option      [-] : valid for layered model only; linear or constant
  @param soil_storage_model_depth          [m] : depth of the soil storage reservoir
  @param soil_moisture_layered_bmi         [-] : if true, soil moisture content of wetting fronts is set by the bmi
  @param soil_moisture_fraction            [-] : fraction of soil moisture in the top 40 cm or user-defined
                                                 soil_moisture_fraction_depth (water in the topsoil over total water)
  @param soil_moisture_fraction_depth      [m] : user specified depth for the fraction of soil moisture (default is 40 cm)

  @param Qb_topmodel           [m/hr] : baseflow in the topmodel
  @param Qv_topmodel           [m/hr] : recharge rate of the saturated zone from the unsaturated zone in the topmodel
  @param global_deficit        [m]    : catchment soil moisture deficit in the topmodel
  @param field_capacity        [m]    : soil field capacity
  @param cat_area              [m^2]  : catchment area
  @param verbosity             [-]    : flag for screen outputs for debugging, options = none, high
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
    int    shape[3];
    double spacing[8];
    double origin[3];

    double soil_storage;
    double soil_storage_change_per_timestep;
    double water_table_depth;
    double *soil_moisture_profile;
    
    double *smcmax;
    double b;
    double satpsi;
    int    ncells;
    double soil_depth;
    double soil_depth_NWM;
    double last_layer_depth;
    double *soil_z;

    double soil_moisture_fraction;
    double soil_moisture_fraction_depth;
    
    int    soil_storage_model;
    double soil_storage_model_depth;
    int    soil_moisture_profile_option;

    bool   init_profile;
    std::string verbosity;
    
    // layered model
    double *soil_moisture_wetting_fronts;
    double *soil_depth_wetting_fronts;
    double *soil_depth_layers;
    int     num_wetting_fronts;
    int     max_num_wetting_fronts;
    int     num_layers;
    bool    soil_depth_layers_bmi;
    bool    smcmax_bmi;

    //topmodel bmi outputs
    double Qb_topmodel;
    double Qv_topmodel;
    double global_deficit;
    double field_capacity;
    int    water_table_based_method;
    double cat_area;
    
  };


  void SoilMoistureProfile(std::string config_file, struct soil_profile_parameters* parameters);

  void InitFromConfigFile(std::string config_file, struct soil_profile_parameters* parameters);

  // reading 1D array from the config file
  std::vector<double> ReadVectorData(std::string param_name, std::string param_value);

  // update the profile for the current timestep
  void SoilMoistureProfileUpdate(struct soil_profile_parameters* parameters);

  // computes soil moisture profile for conceptual reservoir
  void SoilMoistureProfileFromConceptualReservoir(struct soil_profile_parameters* parameters);

  // computes soil moisture profile for layered-reservoir
  void SoilMoistureProfileFromLayeredReservoir(struct soil_profile_parameters* parameters);

  // computes soil moisture profile for Topmodel
  void SoilMoistureProfileFromWaterTableDepth(struct soil_profile_parameters* parameters);
  
  // computes linearly interpolated values for layered-reservoir with option = linear
  double LinearInterpolation(double z, double z1, double z2, double t1, double t2);

  void FindWaterTableLayeredReservoir(struct soil_profile_parameters* parameters);
  // print soil moisture profile and soil depths
  void PrintSoilMoistureProfile(struct soil_profile_parameters* parameters);
  
};

#endif
