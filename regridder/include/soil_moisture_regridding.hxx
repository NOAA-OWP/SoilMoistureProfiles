#ifndef SMCM_H_INCLUDED
#define SMCM_H_INCLUDED

/*
  @author: Ahmad Jan (ahmad.jan@noaa.gov)
  Soil moisture mapping module maps topographic wetness index (TWI) based soil moisture to National Water Model (NWM) 1x1 km grid
  IDEA: Module computes and maps local (subcatchment) soil moisture from catchment to NWM grid (1x1 km).
  The idea of computing local soil moisture is based on the concept of TopModel, where given the global (mean)
  watershed soil moisture deficit and TWI, local soil moisture deficit is computed. Here, we use the watershed
  maximum storage capacity and TWI to compute local soil moisture content which is then mapped onto NWM 1x1 km
  grid using the spatial mapping information provided in the hydrofabric.


  Inputs:
  @param maxsmc  [-]                  : maximum soil moisture (porosity)
  @param cat_global_deficit [m]       : Basin soil moisture deficit
  @param cat_local_deficit [m]        : catchment local soil deficit
  @param cat_local_moisture [m]       : catchment local soil moisture
  @param depth                        : depth/thickness of the soil column (needed to compute maximum soil moisture storage)
  @param szm                          : a decay factor of transmissivity with increase in storage deficit with dimensions of length (m in the TopModel)
  @param grid_id                      : global grid id of the NWM (provided by the hydrofabric)
  @param cat_grid_id                  : catchment id in the spatial data file corresponding to each grid_id. Spatial file has the mapping between catchments' ids and NWM grid
  @param cat_id                       : catchment id in the attribute data file. Length is the number of catchments
  @param grid_id_unique               : unique grid id (removes all duplicate ids from grid_id). 
  @param grid_id_unique_index         : index of the unique grid id. Needed to store the soil moisture fraction at the right index. Note grid_ids are not unique and not sorted and can start anywhere, so we index them based on their first appearance. For example, if grid_id = [19,22,24,31,19,16,9,31] then unique grid id index would be grid_id_unique_index = [0,1,2,4,0,5,6,3]
  @param cat_id_index                 : index of catchment ids. Again since the catchments are not properly ordered so we keep track of them through their indices 
  @param grid_area_fraction  [mxm]    : fractional area of NWM grid contained withing a catchment
  @param TWI                          : topographic wetness index of each catchment
  @param dist_area_TWI                : the distribution of area corresponding to ln(A/tanB) histogram
  @param grid_soil_moisture           : grided soil moisutre content
  @param grid_total_area              : total area of grid. needed for grid other than 1x1
  @param ngrids                       : total number of grid_id (include duplicate ids)
  @param ngrids_unique                : total number of unique grid_id
  @param ncats                        : total number of catchments
  @param cat_storage_max              : catchment maximum soil moisture storage = porosity * depth
  @param cat_global_storage           : basin (all catchments) soil moisture storage
  @param total_area                   :
  @param areal_avg_TWI                : areal (area weighted) average TWI (A_i * TWI_i, i=catchment)
 */

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <numeric>

using namespace std;

namespace smc_mapping {

  class SoilMoistureMapping {
  private:
    //std::string config_file;
        
  public:
    int shape[3];
    int shape_basin[3];
    double spacing[8];
    double origin[3];
    
    double maxsmc; // maximum soil moisture (porosity)
    double cat_global_deficit; // sbar [m]
    double depth; // soil depth/thickness [m]
    double szm; // famous m parameter

    double *cat_local_deficit; // soil moisture deficit [m] (1D array)
    double *cat_local_moisture; // soil moisture storage [m] (1D array)

    /*** 1D arrays *****************/
    int *grid_id;                
    int *cat_grid_id;            //this is the cat ID in the spatial data file (corresponds to each grid id
    int *cat_id;                 // this is the ID in the attribute data file
    int *grid_id_unique;
    int *grid_id_unique_index;
    int *cat_id_unique;
    int *cat_id_unique_index;
    int *cat_id_index; // check if this is needed
    double *grid_area_fraction; 
    double *TWI; 
    double *dist_area_TWI; 
    double *grid_soil_moisture; 
    double *grid_total_area; 
    /******************************/

    /*** 2D arrays ***************/

    double **TWI_V1;
    double **dist_area_TWI_V1;
    double *areal_avg_TWI_V1; // area weighted average TWI
    int num_twi_per_cat;
    double *maxsmc_V1; // maximum soil moisture (porosity)
    double *cat_global_deficit_V1; // sbar [m]
    double *depth_V1; // soil depth/thickness [m]
    double *szm_V1; // famous m parameter
    std::string model_type;
    /******************************/

    int ngrids;
    int length; // total number of data points given in the file
    //    int ngrids_unique; // u: unique
    int ncats;

    double cat_storage_max;
    double cat_global_storage;
    double total_area;
    double areal_avg_TWI; // area weighted average TWI
    
    SoilMoistureMapping();
    SoilMoistureMapping(std::string config_file);

    void SoilMoistureFromBasinToGrid();
    void InitFromConfigFile(std::string config_file);
    void ReadSpatialData(std::string spatial_file);
    void ReadTWIData(std::string spatial_file);
    void ReadCatchmentParamsData(std::string spatial_file);
    void AreaWeightedAverageTWI();
    void ComputeLocalSoilMoisture();
    void ComputeGriddedSoilMoisture();

    ~SoilMoistureMapping();
  };

};


#endif
