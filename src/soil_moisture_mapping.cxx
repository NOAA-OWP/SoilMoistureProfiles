#ifndef SMCMAP_C_INCLUDED
#define SMCMAP_C_INCLUDED

#include <cstring>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include "../include/soil_moisture_mapping.hxx"


/*
Default constructor. 
@param shape, spacing etc. may not be needed in this BMI
*/
smc_mapping::SoilMoistureMapping::
SoilMoistureMapping()
{
  this->shape[0] = 1;
  this->shape[1] = 1;
  this->shape[2] = 1;
  this->shape_basin[0] = 1;
  this->shape_basin[1] = 1;
  this->shape_basin[2] = 1;
  this->spacing[0] = 1.;
  this->spacing[1] = 1.;
  this->origin[0] = 0.;
  this->origin[1] = 0.;
  
}

/*
Constructor called by BMI initialize to read configuration file and initialized parameters
We call local soil moisture and gridded soil moisture methods here to populate them at time 0
*/
smc_mapping::SoilMoistureMapping::
SoilMoistureMapping(std::string config_file)
{
  this->InitFromConfigFile(config_file);
  
  this->shape[0] = 1;
  this->shape[1] = 1;
  this->shape[2] = 1;
  this->shape_basin[0] = 1;
  this->spacing[0] = 1.;
  this->spacing[1] = 1.;
  this->origin[0] = 0.;
  this->origin[1] = 0.;

  this->cat_local_deficit = new double[this->ncats];
  this->cat_local_moisture = new double[this->ncats];

  this->grid_soil_moisture = new double[ngrids_unique];
  this->grid_total_area = new double[ngrids_unique];
  
  this->cat_storage_max = this->depth * this->maxsmc;

  // this needs to be called once; so calling here in the constructor
  AreaWeightedAverageTWI();

  // call to initialize local soil moisture content
  ComputeLocalSoilMoisture();

  // call to map the initialized soil moisture to the NWM grid
  ComputeGridedSoilMoisture();
  
}

/*
BMI call this method to update and map soil moisture from catchment from NWM grid
*/
void smc_mapping::SoilMoistureMapping::
SoilMoistureFromBasinToGrid()
{
  //update global storage
  // this->cat_global_storage = cat_storage_max - *cat_global_deficit;

  //update catchment soil moisture
  ComputeLocalSoilMoisture();

  // update/map grided soil moisture
  ComputeGridedSoilMoisture();
  
}

/*
Computes the areal (area weighted) average of TWI
areal average: 1/Area * integral (sub_cat_area_i * TWI_i)

*/
void smc_mapping::SoilMoistureMapping::
AreaWeightedAverageTWI()
{
  // Get total area (TWI area which should sum to 1.0)
  // and areal average (area weighted-average of TWI) 

  this->areal_avg_TWI = 0.0; 
  
  std::vector<double> dist_area_lnaotb(ncats);
  total_area = 1.0; // hacked value for testing
  
  for (int i=0; i < ncats; i++)
    dist_area_lnaotb[i] = dist_area_TWI[i]/total_area;
  
  for (int i=1; i < ncats; i++)
    this->areal_avg_TWI +=  dist_area_lnaotb[i] * (TWI[i]+TWI[i-1])/2.0;

  this->areal_avg_TWI = 5.454978310883997; //hacked value for testing

}

/*
Computes catchment local deficit and catchment local moisture
@param cat_global_deficit : catchment global deficit provided in the config file
@param szm                : decay factor of transmissivity (hydraulic conductivity * area) with increase in soil deficit
@param areal_avg_TWI      : areal (weighted-average) topographic wetness index
@param TWI                : topographic wetness index
@param cat_local_deficit  : the catchment local deficit (same as computed in TOPMODEL)
@param cat_local_moisture : the catchment local moisture bases on global storage and local deficit
@param ncats              : number of catchments
*/
void smc_mapping::SoilMoistureMapping::
ComputeLocalSoilMoisture()
{
  for (int i=0; i < ncats; i++) {
    cat_local_deficit[i] = cat_global_deficit + szm * (areal_avg_TWI-TWI[i]);
    cat_local_moisture[i] = cat_storage_max - cat_local_deficit[i];
  }
  
}

/*
Computes and maps catchment local moisture to NWM grid
@param ngrids_unique         : total number of unique grid ids
@param ngrids                : total number of grid ids (may include duplicates; in case a grid lies in multiple catchments)
@param grid_soil_moisture    : 1D array of gridded soil moisture
@param grid_total_area       : 1D array of grid total area (if different than 1x1 km) stores gridded soil moisture
@param grid_id_unique_index  : index of the unique grid ids (as the grids ids are neither sorted nor start from 0, we access them bases on their indices)
@param cat_grid_id           : catchment id corresponding to each grid id (this is different than 1D array of catchment ids) 
@param cat_id_index          : index of the unique cat ids (as the catchment ids are neither sorted nor start from 0, we access them bases on their indices)
@param grid_soil_moisture    : grided soil moisutre content
@param grid_area_fraction    : fractional area of NWM grid contained withing a catchment
*/
void smc_mapping::SoilMoistureMapping::
ComputeGridedSoilMoisture()
{
  
  // make sure smc and area are set to zero before mapping
  for (int i=0; i<ngrids_unique;i++) {
    this->grid_soil_moisture[i] = 0.0;
    this->grid_total_area[i] = 0.0;
  }
  
  for (int i=0; i <ngrids; i++) {
    int gid_index = this->grid_id_unique_index[i];
    int map_cid = this->cat_grid_id[i]; //catchment id
    int cat_index = this->cat_id_index[map_cid];

     grid_soil_moisture[gid_index] += grid_area_fraction[i] * cat_local_moisture[cat_index];

     grid_total_area[gid_index] += grid_area_fraction[i];
     
  }

  // std::cout<<"ID,SMC"<<"\n";
  for (int i=0; i <ngrids_unique; i++) {
    //int id = this->grid_id_unique[i];
    grid_soil_moisture[i] /= grid_total_area[i];
    //std::cout<<id<<","<<grid_SMC[i]<<"\n";
  }
  
}

/*
Initialize parameters from the configuration file
*/
void smc_mapping::SoilMoistureMapping::
InitFromConfigFile(std::string config_file)
{ 
  std::ifstream fp;
  fp.open(config_file);
  
  bool is_spatial_data_set = false;
  bool is_TWI_data_set = false;
  bool is_maxsmc_set = false;
  bool is_global_deficit_set = false;
  bool is_depth_set = false;
  bool is_szm_set = false;
  
  while (fp) {
    
    std::string key;
    std::getline(fp, key);
    
    int loc = key.find("=");
    std::string key_sub = key.substr(0,loc);
    
    if (key_sub == "spatial_map_file") {
      std::string tmp_key = key.substr(loc+1,key.length());
      ReadSpatialData(tmp_key);
      is_spatial_data_set = true;
      continue;
    }
    else if (key_sub == "TWI_file") {
      std::string tmp_key = key.substr(loc+1,key.length());
      ReadTWIData(tmp_key);
      is_TWI_data_set = true;
      continue;
    }
    else if (key_sub == "global_deficit") {
      this->cat_global_deficit = std::stod(key.substr(loc+1,key.length()));
      is_global_deficit_set = true;
      continue;
    }
    else if (key_sub == "depth") {
      this->depth = std::stod(key.substr(loc+1,key.length()));
      assert (this->depth > 0);
      is_depth_set = true;
      continue;
    }
    else if (key_sub == "porosity") {
      this->maxsmc = std::stod(key.substr(loc+1,key.length()));
      is_maxsmc_set = true;
      continue;
    }
    else if (key_sub == "szm") {
      this->szm = std::stod(key.substr(loc+1,key.length()));
      is_szm_set = true;
      continue;
    }
    
  }
  fp.close();
  
  if (!is_spatial_data_set) {
    std::stringstream errMsg;
    errMsg << "Spatail data file not provided in the config file "<< config_file << "\n";
    throw std::runtime_error(errMsg.str());
  }

  if (!is_maxsmc_set) {
    std::stringstream errMsg;
    errMsg << "Porosity not set in the config file "<< config_file << "\n";
    throw std::runtime_error(errMsg.str());
  }
  
  if (!is_global_deficit_set) {
    std::stringstream errMsg;
    errMsg << "Global deficit not set in the config file "<< config_file << "\n";
    throw std::runtime_error(errMsg.str());
  }
  
  if (!is_szm_set) {
    std::stringstream errMsg;
    errMsg << "Parameter szm not set in the config file "<< config_file << "\n";
    throw std::runtime_error(errMsg.str());
  }

}

/*
Reads the spatial file that provides mapping between grid ids (NWM) and catchment ids and the area fraction of each grid within the catchment

*/
void smc_mapping::SoilMoistureMapping::
ReadSpatialData(std::string spatial_file)
{
  std::ifstream fp;
  fp.open(spatial_file);
  if (!fp) {
    cout<<"file "<<spatial_file<<" doesn't exist. \n";
    abort();
  }

  // vectors with _v are local vectors created to perform basic operations that we can't do with arrays directly (e.g., sorting)
  std::vector<int> grid_id_v(0.0);
  std::vector<int> cat_id_v(0.0);
  std::vector<double> area_fraction_v(0.0);
  std::vector<double> TWI_v(0.0);
  std::vector<double> dist_area_TWI_v(0.0);
  std::vector<double> local_storage_v(0.0);
  std::vector<string> vars;
  std::string line, cell;
  
  //read first line of strings which contains variables names.
  std::getline(fp, line);
  std::stringstream lineStream(line);
 
  while(std::getline(lineStream,cell, ',')) {
    vars.push_back(cell);
  }

  // check the order of the variables; this can be avoided
  if (vars[0].compare("grid_id") != 0) {
    std::stringstream errMsg;
    errMsg << "data order issue: grid_id should be the first column "<<"\n";
    throw std::runtime_error(errMsg.str());
  }
  else if (vars[1].compare("cat_id") != 0) {
    std::stringstream errMsg;
    errMsg << "data order issue: cat_id should be the second column "<<"\n";
    throw std::runtime_error(errMsg.str());
  }
  else if (vars[2].compare("area_fraction") != 0) {
    std::stringstream errMsg;
    errMsg << "data order issue: area_fraction should be the third column "<<"\n";
    throw std::runtime_error(errMsg.str());
  }

  int count =0;
  while (fp) {
    std::getline(fp, line);
    std::stringstream lineStream(line);
    
    while(std::getline(lineStream,cell, ',')) {
      
      if (count % 3 == 0)
	grid_id_v.push_back(stod(cell));
      if (count % 3 == 1)
	cat_id_v.push_back(stod(cell));
      if (count % 3 == 2)
	area_fraction_v.push_back(stod(cell));
      count++;
    }
    
  }
  
  int size_v = grid_id_v.size();

  // memory allocation
  this->grid_id = new int[size_v];
  this->cat_grid_id = new int[size_v];
  this->grid_area_fraction = new double[size_v];

  // copyback data from vectors to 1D arrays
  for (int i=0; i<size_v; i++) {
    this->grid_id[i] = grid_id_v[i];
    this->cat_grid_id[i] = cat_id_v[i] - 1;
    this->grid_area_fraction[i] = area_fraction_v[i];
  }
  
  this->ngrids = size_v;
  
  //std::vector<int> cat_id_tmp = cat_id_v;
  //std::sort(cat_id_tmp.begin(), cat_id_tmp.end());
  //*this->ncats = std::unique(cat_id_tmp.begin(), cat_id_tmp.end()) - cat_id_tmp.begin();


  // We need unique grid ids for storing soil moisture data 
  std::vector<int> grid_id_index_v = grid_id_v;
  std::vector<int>::iterator begin = grid_id_index_v.begin();
  std::vector<int>::iterator end = grid_id_index_v.end();

  // remove duplicate ids
  for (std::vector<int>::iterator it = begin; it != end; ++it) {
    end = std::remove(it + 1, end, *it); 
  }
  
  grid_id_index_v.erase(end, grid_id_index_v.end()); // remove does not actually remove it, but put it toward the end of the array, so let's erase all those entities

  this->ngrids_unique = grid_id_index_v.size();
  
  this->grid_id_unique = new int[ngrids_unique];

  for (int i=0; i<ngrids_unique; i++) {
    this->grid_id_unique[i] = grid_id_index_v[i];
  }

  // a temporary vector of size ngrids to keep unique index of each grid_id
  // grid_id = [19,22,24,31,19,16,9,31] unique index would be grid_id_unique_index = [0,1,2,4,0,5,6,3]
  std::vector<int> vec(ngrids);

  // Get index for unique ids
  for (int i=0; i<ngrids_unique; i++) { //start with unique id
    for (int j=0; j<ngrids; j++) { // loop over NWM grid id to cover duplicate ids and assign them same index
      if (grid_id_unique[i] == grid_id[j]) {
	vec[j] = i;
      }
    }
  }
  
  this->grid_id_unique_index = new int[ngrids];

  // copy vector entities into 1D array
  for (int i=0; i<ngrids; i++) {
    this->grid_id_unique_index[i] = vec[i];
  }
  
}

/*
Reads catchment ids, topographic wetness index, and distribution of area corresponding to TWI histrogram
*/
void smc_mapping::SoilMoistureMapping::
ReadTWIData(std::string spatial_file)
{
  std::ifstream fp;
  fp.open(spatial_file);
  if (!fp) {
    std::stringstream errMsg;
    errMsg << "file "<<spatial_file<<" doesn't exist. "<<"\n";
    throw std::runtime_error(errMsg.str());
  }

  std::vector<int> cat_id_v(0.0);
  std::vector<double> TWI_v(0.0);
  std::vector<double> dist_area_TWI_v(0.0);
  std::vector<string> vars;
  std::string line, cell;
  
  //read first line of strings that contains variables names.
  std::getline(fp, line);
  std::stringstream lineStream(line);
    
  while(std::getline(lineStream,cell, ',')) {
    vars.push_back(cell);
  }

  
  if (vars[0].compare("cat_id") != 0) {
    std::stringstream errMsg;
    errMsg << "data order issue: cat_id should be the 1st column "<<"\n";
    throw std::runtime_error(errMsg.str());
  }
  else if (vars[1].compare("TWI") != 0) {
    std::stringstream errMsg;
    errMsg << "data order issue: TWI should be the 2nd column "<<"\n";
    throw std::runtime_error(errMsg.str());
  }
  else if (vars[2].compare("dist_area_TWI") != 0) {
    std::stringstream errMsg;
    errMsg << "data order issue: distributed area TWI should be the 3rd column "<<"\n";
    throw std::runtime_error(errMsg.str());
  }

  int count =0; 
  while (fp) {
    std::getline(fp, line);
    std::stringstream lineStream(line);
    
    while(std::getline(lineStream,cell, ',')) {
      
      if (count % 3 == 0)
	cat_id_v.push_back(stod(cell));
      if (count % 3 == 1)
	TWI_v.push_back(stod(cell));
      if (count % 3 == 2)
	dist_area_TWI_v.push_back(stod(cell));
      count ++;
    }
    
  }
  
  int size_v = cat_id_v.size();
 
  this->cat_id = new int[size_v];
  this->TWI = new double[size_v];
  this->dist_area_TWI = new double[size_v];

  // copy vector data to 1D array (1D arrays are data members)
  for (int i=0; i<size_v; i++) {
    this->cat_id[i] = cat_id_v[i]; // check if cat_id start at 0 or 1.
    this->TWI[i] = TWI_v[i];
    this->dist_area_TWI[i] = dist_area_TWI_v[i];
  }


  // total area of the distribution of areal TWI; it should sum to 1.0, computing it here would be good to checking 
  total_area = std::accumulate(dist_area_TWI_v.begin(), dist_area_TWI_v.end(), 0.0);

  this->cat_id_index = new int[size_v];

  // correct index is need to extract the right catchment id as the catchment ids may not be in order
  this->ncats = size_v;
  for (int i=0; i<ncats; i++) { 
    this->cat_id_index[cat_id[i]-1] = i;   
  }
  
}

smc_mapping::SoilMoistureMapping::
~SoilMoistureMapping()
{}

#endif
