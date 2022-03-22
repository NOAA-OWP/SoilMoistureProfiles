#ifndef SMCM_H_INCLUDED
#define SMCM_H_INCLUDED


#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <numeric>

using namespace std;

namespace smc_mapping {

  class SMCMapping {
  private:
    //std::string config_file;
        
  public:
    int shape[3];
    int shape_basin[3];
    double spacing[8];
    double origin[3];
    
    //    double *soil_texture;
    double *phi; // porosity
    double *cat_global_deficit; // sbar [m]
    double *cat_local_deficit; // soil moisture deficit [m]
    double *cat_local_moisture; // soil moisture storage [m]
    double *depth; // soil depth/thickness [m]
    double *szm; // famous m parameter
    
    
    int *grid_id;
    int *cat_grid_id; //this is the cat ID in the spatial data file (corresponds to each grid id)
    int *cat_id; // this is the ID in the attribute data file
    int *grid_id_unique;
    int *grid_id_unique_index;
    int *cat_id_index;
    double *grid_area_fraction;
    double *TWI;
    double *dist_area_TWI;
    double *grid_SMC;
    double *grid_total_area;
    
    int *ngrids;
    int *ngrids_u; // u: unique
    int *ncats;

    double cat_storage_max;
    double cat_global_storage;
    double total_area;
    double areal_avg_TWI; // area weighted average TWI
    
    SMCMapping();
    SMCMapping(std::string config_file);

    void SMCFromBasinToGrid();
    void InitFromConfigFile(std::string config_file);
    void ReadSpatialData(std::string spatial_file);
    void ReadTWIData(std::string spatial_file);
    void AreaWeightedAverageTWI();
    void ComputeLocalSoilMoisture();
    void ComputeGridedSoilMoisture();
    ~SMCMapping();
  };

};


#endif
