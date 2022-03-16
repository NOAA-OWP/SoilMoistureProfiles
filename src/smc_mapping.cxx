#ifndef SMCMAP_C_INCLUDED
#define SMCMAP_C_INCLUDED

#include <cstring>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include "../include/smc_mapping.hxx"

smc_mapping::SMCMapping::
SMCMapping()
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

smc_mapping::SMCMapping::
SMCMapping(std::string config_file)
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

  this->cat_local_deficit = new double[*this->ncats];
  this->cat_local_moisture = new double[*this->ncats];

  this->grid_SMC = new double[*ngrids_u];

  this->cat_storage_max = (*this->depth) * (*this->phi);
  
  ComputeLocalMoisture();
  ComputeGridedMoisture();
}

void smc_mapping::SMCMapping::
ComputeLocalMoisture()
{

  // cat_global_deficit updates every timestep
  this->cat_global_storage = cat_storage_max - *cat_global_deficit;

  // Get total area (TWI area which should sum to 1.0)
  // and areal average (area weighted-average of TWI) 

  double TL = 0.0; // areal average: 1/Area * integral (sub_cat_area_i * TWI_i)

  std::vector<double> dist_area_lnaotb(*ncats);
  total_area = 1.0; // hacked value for testing
  
  for (int i=0; i < *ncats; i++)
    dist_area_lnaotb[i] = dist_area_TWI[i]/total_area;

  for (int i=1; i < *ncats; i++)
    TL +=  dist_area_lnaotb[i]*(TWI[i]+TWI[i-1])/2.0;

  TL = 5.454978310883997; //hacked value for testing

  for (int i=0; i < *ncats; i++) {
    cat_local_deficit[i] = (*cat_global_deficit) + (*szm) * (TL-TWI[i]);
    cat_local_moisture[i] = cat_storage_max - cat_local_deficit[i];

  }
  
}

void smc_mapping::SMCMapping::
ComputeGridedMoisture()
{
  double *grid_total_area = new double[*ngrids_u];
  
  for (int i=0; i <*ngrids; i++) {
    int gid_index = this->grid_id_unique_index[i];
    int map_cid = this->cat_grid_id[i]; //catchment id
    int cat_index = this->cat_id_index[map_cid];

     grid_SMC[gid_index] = grid_SMC[gid_index] + grid_area_fraction[i] * cat_local_moisture[cat_index];

     grid_total_area[gid_index] = grid_total_area[gid_index] + grid_area_fraction[i];
     
  }

  std::cout<<"ID,SMC"<<"\n";
  for (int i=0; i <*ngrids_u; i++) {
    int id = this->grid_id_unique[i];
    grid_SMC[i] /= grid_total_area[i];
    std::cout<<id<<","<<grid_SMC[i]<<"\n";
  }
}


void smc_mapping::SMCMapping::
InitFromConfigFile(std::string config_file)
{ 
  std::ifstream fp;
  fp.open(config_file);
  
  bool is_spatial_data_set = false;
  bool is_TWI_data_set = false;
  bool is_phi_set = false;
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
      this->cat_global_deficit = new double[1];
      *this->cat_global_deficit = std::stod(key.substr(loc+1,key.length()));
      is_global_deficit_set = true;
      continue;
    }
    else if (key_sub == "depth") {
      this->depth = new double[1];
      *this->depth = std::stod(key.substr(loc+1,key.length()));
      assert (*this->depth > 0);
      is_depth_set = true;
      continue;
    }
    else if (key_sub == "porosity") {
      this->phi = new double[1];
      *this->phi = std::stod(key.substr(loc+1,key.length()));
      is_phi_set = true;
      continue;
    }
    else if (key_sub == "szm") {
      this->szm = new double[1];
      *this->szm = std::stod(key.substr(loc+1,key.length()));
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

  if (!is_phi_set) {
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

void smc_mapping::SMCMapping::
ReadSpatialData(std::string spatial_file)
{
  std::ifstream fp;
  fp.open(spatial_file);
  if (!fp) {
    cout<<"file "<<spatial_file<<" doesn't exist. \n";
    abort();
  }

  std::vector<int> grid_id_v(0.0);
  std::vector<int> cat_id_v(0.0);
  std::vector<double> area_fraction_v(0.0);
  std::vector<double> TWI_v(0.0);
  std::vector<double> dist_area_TWI_v(0.0);
  std::vector<double> local_storage_v(0.0);
  std::vector<string> vars;
  std::string line, cell;
  
  //read first line of strings which contains forcing variables names.
  std::getline(fp, line);
  std::stringstream lineStream(line);
    
  while(std::getline(lineStream,cell, ',')) {
    vars.push_back(cell);
  }


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
 
  this->grid_id = new int[size_v];
  this->cat_grid_id = new int[size_v];
  this->grid_area_fraction = new double[size_v];
  this->ngrids = new int[1];
  this->ngrids_u = new int[1];
  
  for (int i=0; i<size_v; i++) {
    this->grid_id[i] = grid_id_v[i];
    this->cat_grid_id[i] = cat_id_v[i] - 1;
    this->grid_area_fraction[i] = area_fraction_v[i];
  }
  
  *this->ngrids = size_v;
  
  //std::vector<int> cat_id_tmp = cat_id_v;
  //std::sort(cat_id_tmp.begin(), cat_id_tmp.end());
  //*this->ncats = std::unique(cat_id_tmp.begin(), cat_id_tmp.end()) - cat_id_tmp.begin();


  // We will need unique grid ids for storing soil moisture data 
  std::vector<int> grid_id_index_v = grid_id_v;
  std::vector<int>::iterator begin = grid_id_index_v.begin();
  std::vector<int>::iterator end = grid_id_index_v.end();
  
  for (std::vector<int>::iterator it = begin; it != end; ++it) {
    end = std::remove(it + 1, end, *it);
  }
  
  grid_id_index_v.erase(end, grid_id_index_v.end());

  *this->ngrids_u = grid_id_index_v.size();
  
  this->grid_id_unique = new int[*ngrids_u];

  for (int i=0; i<*ngrids_u; i++) {
    this->grid_id_unique[i] = grid_id_index_v[i];
  }


  this->grid_id_unique_index = new int[*ngrids];
  //std::vector<int> vec(*ngrids);
  
  for (int i=0; i<*ngrids_u; i++) { //start with unique id
    for (int j=0; j<*ngrids; j++) { // loop over nwm grid to cover duplicate ids and assign them same index
      if (grid_id_unique[i] == grid_id[j]) {
	//vec[j] = i;
	this->grid_id_unique_index[i] = vec[i];
      }
    }
  }
  
  /*
  for (int i=0; i<*ngrids; i++) {
    this->grid_id_unique_index[i] = vec[i];
    }*/

}


void smc_mapping::SMCMapping::
ReadTWIData(std::string spatial_file)
{
  std::ifstream fp;
  fp.open(spatial_file);
  if (!fp) {
    cout<<"file "<<spatial_file<<" doesn't exist. \n";
    abort();
  }

  std::vector<int> cat_id_v(0.0);
  std::vector<double> TWI_v(0.0);
  std::vector<double> dist_area_TWI_v(0.0);
  std::vector<string> vars;
  std::string line, cell;
  
  //read first line of strings which contains forcing variables names.
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
  this->ncats = new int[1];
  
  for (int i=0; i<size_v; i++) {
    this->cat_id[i] = cat_id_v[i]; // check if cat_id start at 0 or 1.
    this->TWI[i] = TWI_v[i];
    this->dist_area_TWI[i] = dist_area_TWI_v[i];
  }

  
  total_area = std::accumulate(dist_area_TWI_v.begin(), dist_area_TWI_v.end(), 0.0);

  this->cat_id_index = new int[size_v];

  // correct index is need to extract the right catchment id as the catchments may not be in order
  *this->ncats = size_v;
  for (int i=0; i<*ncats; i++) { 
    this->cat_id_index[cat_id[i]-1] = i;   
  }
  
}

void smc_mapping::SMCMapping::
MapFromBasinToGrid()
{}

smc_mapping::SMCMapping::
~SMCMapping()
{}

#endif
