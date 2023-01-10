#ifndef SMCMAP_C_INCLUDED
#define SMCMAP_C_INCLUDED

#include <cstring>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include "../include/soil_moisture_regridding.hxx"


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

  std::cout<<"Model type : "<<model_type<<"\n";
  std::cout<<"ncats, ngrids "<<ncats<<" "<<ngrids<<"\n";
  
  this->grid_soil_moisture = new double[ngrids]; // total number to NWM grids
  this->grid_total_area = new double[ngrids];
  this->cat_local_moisture = new double[ncats]; // total number of catchments

  
  if (model_type == "soil_deficit_based") {
    this->cat_local_deficit = new double[this->ncats];
      
    this->cat_storage_max = this->depth * this->maxsmc;
    
    // this needs to be called once; so calling here in the constructor
    AreaWeightedAverageTWI();
    
    // call to initialize local soil moisture content
    ComputeLocalSoilMoisture();
    
    // call to map the initialized soil moisture to the NWM grid
    ComputeGriddedSoilMoisture();
  }
  else if (model_type == "soil_moisture_based") {
    // synthetic test case2
    //double syn_smc[] = {0.29,0.03,0.51,0.7,0.09,0.66,0.37,0.44,0.1,0.85,0.7,0.61,0.76,0.48,0.75,0.44,0.84};
    // synthetic test case3
    double syn_smc[] = {0.29,0.03,0.51,0.7,0.09,0.66,0.37,0.44,0.1,0.85,0.7,0.61, 0.48,0.75,0.84, 0.76, 0.44}; 
    
    for (int i=0; i<ncats; i++)
      cat_local_moisture[i] = syn_smc[i];
    
    ComputeGriddedSoilMoisture();
  }
  
}

/*
BMI call this method to update and map soil moisture from catchment from NWM grid
*/
void smc_mapping::SoilMoistureMapping::
SoilMoistureFromBasinToGrid()
{
  //update global storage
  // this->cat_global_storage = cat_storage_max - *cat_global_deficit;

  if (model_type == "soil_deficit_based") {
    //update catchment soil moisture
    ComputeLocalSoilMoisture();

  }
  else if (model_type == "soil_moisture_based") {
    // update/map grided soil moisture
    ComputeGriddedSoilMoisture();
  }
  
}

/*
  For each catchment, compute the areal (area weighted) average of TWI
  areal average: 1/Area * integral (sub_cat_area_i * TWI_i)
  Note: the way hydrofabric generates frequencies (i.e., dist_area_twi) they sum to 1, so Area = 1.
*/

void smc_mapping::SoilMoistureMapping::
AreaWeightedAverageTWI()
{
  // Get total area (TWI area which should sum to 1.0)
  // and areal average (area weighted-average of TWI) 
  /*
  this->areal_avg_TWI = 0.0;

  
  std::vector<double> dist_area_lnaotb(ncats);
  total_area = 1.0; // hacked value for testing
  
  for (int i=0; i < ncats; i++)
    dist_area_lnaotb[i] = dist_area_TWI[i]/total_area;
  
  for (int i=1; i < ncats; i++)
    this->areal_avg_TWI +=  dist_area_lnaotb[i] * (TWI[i]+TWI[i-1])/2.0;

  this->areal_avg_TWI = 5.454978310883997; //hacked value for testing
  */

  this->areal_avg_TWI_V1 = new double[ncats];
  std::vector<double> dist_area_lnaotb(ncats);
  
  double sum_v1;
  for (int i=0; i < ncats; i++) {
    sum_v1 = 0.0;
    for (int j=1; j < num_twi_per_cat; j++) {
      sum_v1 += this->dist_area_TWI_V1[i][j] * (TWI_V1[i][j] + TWI_V1[i][j-1])/2.0;
      //std::cout<<"valu = "<<this->dist_area_TWI_V1[i][j] <<" "<<(TWI_V1[i][j] + TWI_V1[i][j-1])/2.0<<"\n";
      //sum_v1 + = this->dist_area_TWI_V1[i][j] * TWI[i]; // or the simplist approach is to dist_area_twi[i] * twi[i]; also make j=0
    }
    
    //sum1 /= total_area; // should divide by the total dist_area_twi if it does not sum to ONE.
    this->areal_avg_TWI_V1[i] = sum_v1;
   
  }
  std::cerr<<"area_avg = "<<this->areal_avg_TWI_V1[0]<<"\n";
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
  /*
  for (int i=0; i < ncats; i++) {
    cat_local_deficit[i] = cat_global_deficit + szm * (areal_avg_TWI-TWI[i]);
    cat_local_moisture[i] = cat_storage_max - cat_local_deficit[i];
  }
  */
  for (int i=0; i < ncats; i++) {
    for (int j=1; j < num_twi_per_cat; j++) {
    cat_local_deficit[i] = cat_global_deficit_V1[i] + szm_V1[i] * (areal_avg_TWI_V1[i] - TWI_V1[i][j]);
    cat_local_moisture[i] = cat_storage_max - cat_local_deficit[i];
  }

  }
}

/*
Computes and maps catchment local moisture to NWM grid
@param ngrids         : total number of unique grid ids
@param ngrids                : total number of grid ids (may include duplicates; in case a grid lies in multiple catchments)
@param grid_soil_moisture    : 1D array of gridded soil moisture
@param grid_total_area       : 1D array of grid total area (if different than 1x1 km) stores gridded soil moisture
@param grid_id_unique_index  : index of the unique grid ids (as the grids ids are neither sorted nor start from 0, we access them bases on their indices)
@param cat_grid_id           : catchment id corresponding to each grid id (this is different from 1D array of catchment ids) 
@param cat_id_index          : index of the unique cat ids (as the catchment ids are neither sorted nor start from 0, we access them based on their indices)
@param grid_soil_moisture    : grided soil moisutre content
@param grid_area_fraction    : fractional area of NWM grid contained withing a catchment
*/
void smc_mapping::SoilMoistureMapping::
ComputeGriddedSoilMoisture()
{
  std::cout<<"length, ngrids = "<<length<<" "<<ngrids<<"\n";
  
  // make sure smc and area are set to zero before mapping
  for (int i=0; i<ngrids;i++) {
    this->grid_soil_moisture[i] = 0.0;
    this->grid_total_area[i] = 0.0;
  }

  int gid_index, cat_index;
  
  for (int i=0; i <length; i++) {
    
    gid_index = this->grid_id_unique_index[i];
    cat_index = this->cat_id_unique_index[i];
    
    grid_soil_moisture[gid_index] += grid_area_fraction[i] * cat_local_moisture[cat_index];
    
    grid_total_area[gid_index] += grid_area_fraction[i];
    
    //std::cout<<"A2 = "<<i<<" "<<gid_index<<" "<<cat_index<<" "<<cat_local_moisture[cat_index]<<" "<<grid_area_fraction[i]<<"\n";
  }

  std::cout<<"ID,SMC"<<"\n";
  for (int i=0; i <ngrids; i++) {
    int id = this->grid_id_unique[i];
    grid_soil_moisture[i] /= grid_total_area[i];
    std::cout<<id<<","<<grid_soil_moisture[i]<<"\n";
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
  bool is_cat_params_file_set = false;
  bool is_model_type_set = false;
  
  while (fp) {
    
    std::string line;
    string param_key, param_value, param_unit;
    std::getline(fp, line);
    
    int loc_eq = line.find("=") + 1;
    int loc_u = line.find("[");
    param_key = line.substr(0,line.find("="));

    bool is_unit = line.find("[") != string::npos;

    if (is_unit)
      param_unit = line.substr(loc_u,line.find("]")+1);
    else
      param_unit = "";

    param_value = line.substr(loc_eq,loc_u - loc_eq);
    
    if (param_key == "model_type") {
      if (param_value == "soil_moisture_based" || param_value == "Soil_moisture_based")
	model_type = "soil_moisture_based";
      else if (param_value == "soil_deficit_based" || param_value == "Soil_deficit_based")
	model_type = "soil_deficit_based";

      is_model_type_set = true;
      continue;
    }
    else if (param_key == "grid_cat_ids_map_file") {
      ReadSpatialData(param_value);
      is_spatial_data_set = true;
      continue;
    }
    else if (param_key == "TWI_file") {
      ReadTWIData(param_value);
      is_TWI_data_set = true;
      continue;
    }
    else if (param_key == "cat_params_file") {
      ReadCatchmentParamsData(param_value);
      is_cat_params_file_set = true;
      continue;
    }
    else if (param_key == "global_deficit") {
      this->cat_global_deficit = std::stod(param_value); //std::stod(key.substr(loc+1,key.length()));
      is_global_deficit_set = true;
      continue;
    }
    else if (param_key == "depth") {
      this->depth = std::stod(param_value); // std::stod(key.substr(loc+1,key.length()));
      assert (this->depth > 0);
      is_depth_set = true;
      continue;
    }
    else if (param_key == "porosity") {
      this->maxsmc = std::stod(param_value); //std::stod(key.substr(loc+1,key.length()));
      is_maxsmc_set = true;
      continue;
    }
    else if (param_key == "szm") {
      this->szm = std::stod(param_value); //std::stod(key.substr(loc+1,key.length()));
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

  if (!is_model_type_set) {
    std::stringstream errMsg;
    errMsg << "Model type not specified in the config file "<< config_file << "\n";
    throw std::runtime_error(errMsg.str());
  }
  
  if (model_type == "soil_deficit_based") {

    if (!is_TWI_data_set) {
	std::stringstream errMsg;
	errMsg << "TWI data file not provided in the config file "<< config_file << "\n";
      throw std::runtime_error(errMsg.str());
      }
    
    if (!is_cat_params_file_set) {
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

    
    
  }
 
}

/*
  Reads the spatial file that provides mapping between grid ids (NWM) and catchment ids and the area fraction of each grid within the catchment
  the module mainly provides four 1D arrays
  grid_id_unique       : 1D array of unique NWM grid ids, has a size of ngrids
  cat_id_unique        : 1D array of unique catchments ids, has a size of ncats
  grid_id_unique_index : 1D array of unique indices NWM grid ids, has a size of length (total number of data points given)  
  cat_id_unique_index  : 1D array of unique indices catchment ids, has a size of length (total number of data points given)
*/
void smc_mapping::SoilMoistureMapping::
ReadSpatialData(std::string spatial_file)
{
  std::ifstream fp;
  fp.open(spatial_file);

  if (!fp) {
    cout<<"file "<<spatial_file<<" doesn't exist. \n";
    exit(0);
  }

  // vectors with _v are local vectors created to perform basic operations that we can't do with arrays directly (e.g., sorting)
  std::vector<int> grid_id_v(0.0);
  std::vector<int> cat_id_v(0.0);
  std::vector<double> area_fraction_v(0.0);
  std::vector<string> vars;
  std::string line, cell;
  int vars_size = -1;
  
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

  vars_size = vars.size();
  
  int count = 0;
  while (std::getline(fp, line)) {
    std::stringstream lineStream(line);
    
    while(std::getline(lineStream,cell, ',')) {
      
      if (count % vars_size == 0)
	grid_id_v.push_back(stod(cell));
      else if (count % vars_size == 1)
	cat_id_v.push_back(stod(cell));
      else if (count % vars_size == 2)
	area_fraction_v.push_back(stod(cell));

      count++;
      
    }
    
  }
  
  int size_v = grid_id_v.size();
  this->length = size_v; // not unique; may include/count duplicate ids; the total number of data points given in the file
  
  // memory allocation
  this->grid_id = new int[length];
  this->cat_grid_id = new int[length];
  this->grid_area_fraction = new double[length];

  // copyback data from vectors to 1D arrays
  for (int i=0; i<length; i++) {
    this->grid_id[i] = grid_id_v[i];
    this->cat_grid_id[i] = cat_id_v[i];
    this->grid_area_fraction[i] = area_fraction_v[i];
  }
  

  //***********************************************************************************
  // sort catchment ids vector and remove duplicate ids
  std::vector<int> cat_id_tmp = cat_id_v;
  std::sort(cat_id_tmp.begin(), cat_id_tmp.end());

  std::vector<int> cat_id_unique_v = cat_id_v;
  std::vector<int>::iterator begin = cat_id_unique_v.begin();
  std::vector<int>::iterator end   = cat_id_unique_v.end();

  // removing duplicate ids
  for (std::vector<int>::iterator it = begin; it != end; ++it) {
    end = std::remove(it + 1, end, *it); 
  }
  
  cat_id_unique_v.erase(end, cat_id_unique_v.end()); /* remove does not actually remove it, but put it toward the
						      end of the array, so let's erase all those entities */
  this->ncats = cat_id_unique_v.size();

  this->cat_id_unique = new int[ncats];

  // save unique catchment ids to 1D array
  for (int i=0; i<ncats; i++) {
    this->cat_id_unique[i] = cat_id_unique_v[i];
    std::cout<<"cat ids unique = "<<this->cat_id_unique[i]<<"\n";
  }
   
  std::vector<int> vec_v(length);

  // Get index for unique ids
  for (int i=0; i<length; i++) { //loop over length of data
    for (int j=0; j<ncats; j++) { // loop over cat ids to cover duplicate ids and assign them same index
      
      if (cat_id_unique[j] == cat_grid_id[i]) {
	//std::cout<<"Test = "<< i <<" "<<j <<" "<<cat_id_unique[j] <<" "<<cat_grid_id[i]<<"\n";
	vec_v[i] = j;
	//	break;
      }
    }
  }

  // each catchment has a unique index
  this->cat_id_unique_index = new int[length];

  // copy vector entities into 1D array
  // if cat_id = [19,22,24,31,19,16,9,31] then unique index would be cat_id_unique_index = [0,1,2,3,0,4,5,3]
  for (int i=0; i<length; i++) {
    this->cat_id_unique_index[i] = vec_v[i];
    std::cout<<"cat ids unique index = "<<vec_v[i]<<"\n";
  }
  
  //***********************************************************************************
  // repeating the above process for grid ids
  // We need unique grid ids for storing soil moisture data
  
  std::vector<int> grid_id_unique_v = grid_id_v;
  begin = grid_id_unique_v.begin();
  end = grid_id_unique_v.end();

  // remove duplicate ids
  for (std::vector<int>::iterator it = begin; it != end; ++it) {
    end = std::remove(it + 1, end, *it); 
  }
  
  grid_id_unique_v.erase(end, grid_id_unique_v.end());

  this->ngrids = grid_id_unique_v.size();
  
  this->grid_id_unique = new int[ngrids];

  for (int i=0; i<ngrids; i++) {
    this->grid_id_unique[i] = grid_id_unique_v[i];
    std::cout<<"grid_id_un = "<<this->grid_id_unique[i]<<"\n";
  }

  // a temporary vector of size length to keep unique index of each grid_id
  // grid_id = [19,22,24,31,19,16,9,31] unique index would be grid_id_unique_index = [0,1,2,3,0,4,5,3]
  std::vector<int> vec(length);

  // Get index for unique ids
  for (int i=0; i<ngrids; i++) { //start with unique id
    for (int j=0; j<length; j++) { // loop over NWM grid id to cover duplicate ids and assign them same index
      if (grid_id_unique[i] == grid_id[j]) {
	vec[j] = i;
      }
    }
  }
  
  this->grid_id_unique_index = new int[length];

  // copy vector entities into 1D array
  for (int i=0; i<length; i++) {
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

  int count = 0;
  std::getline(fp, line);
  this->num_twi_per_cat = stod(line); // number of bins in TWI distribution (or TWI values)
  //std::cout<<"TWIv1 = "<<line<<"\n";

  std::vector<std::vector<double> > TWI_v1;
  std::vector<std::vector<double> > dist_area_TWI_v1;

  std::string delimiter = ":";
  int loc = -1;
  while (std::getline(fp, line)) {
    
    
    loc = line.find(delimiter);
    //std::cout<<"line = "<<line.substr(0,loc)<<" "<<line.substr(loc+1,std::string::npos) <<"\n";

    
    std::vector<double> v1, v2;
    std::stringstream lineStream(line.substr(loc+1,std::string::npos));

    //v1.push_back(stod(line.substr(0,loc)));
    //v2.push_back(stod(line.substr(0,loc)));
    cat_id_v.push_back(stod(line.substr(0,loc)));
    count++;
    //std::cerr<<"size = "<< cat_id_v.size()<<" "<<loc<<"\n";
    //abort();
    int i = 0;
    while(std::getline(lineStream,cell, ',')) {
      //std::cout<<"val1 = "<< cell<<"\n";
       if (i < num_twi_per_cat)
	 v1.push_back(stod(cell));
       else
	 v2.push_back(stod(cell));
       i++;
    }
    
    dist_area_TWI_v1.push_back(v1);
    TWI_v1.push_back(v2);

  }
  for (int j=0; j< num_twi_per_cat; j++)
    std::cout<<"val = "<<TWI_v1[0][j]<<" "<<dist_area_TWI_v1[0][j]<<"\n";
  //abort();
  /*
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
    std::cout<<"TWIv = "<<cell<<" \n";
    abort();
  }
  */
  std::cout<<"size = "<< cat_id_v.size()<<" "<<count<<"\n";
  int size_v = cat_id_v.size();
 
  this->cat_id = new int[size_v];
  this->TWI = new double[size_v];
  this->dist_area_TWI = new double[size_v];
  
  // copy vector data to 1D array (1D arrays are data members)
  for (int i=0; i<size_v; i++) {
    this->cat_id[i] = cat_id_v[i]; // check if cat_id start at 0 or 1.
    //this->TWI[i] = TWI_v[i];
    //this->dist_area_TWI[i] = dist_area_TWI_v[i];
  }
  
  this->TWI_V1 = new double*[size_v];
  this->dist_area_TWI_V1 = new double*[size_v];
  
  for (int i=0; i<size_v; i++) {
    TWI_V1[i] = new double[num_twi_per_cat];
    dist_area_TWI_V1[i] = new double[num_twi_per_cat];
  }

  for (int i=0; i<size_v; i++) {
    for (int j=0; j<num_twi_per_cat; j++) {
      this->TWI_V1[i][j] = TWI_v1[i][j];
      this->dist_area_TWI_V1[i][j] = dist_area_TWI_v1[i][j];
      //std::cout<<"A1 = "<< this->TWI_V1[0][j] <<" "<<TWI_v1[i][j]<<"\n";
    }
    //    abort();
  }

  // for (int j=0; j<num_twi_per_cat; j++)
  //  std::cout<<"A1 = "<< this->TWI_V1[0][j] <<"\n";
  
  // total area of the distribution of areal TWI; it should sum to 1.0, computing it here would be good for checking 
  total_area = std::accumulate(dist_area_TWI_v.begin(), dist_area_TWI_v.end(), 0.0);

  this->cat_id_index = new int[size_v];

  // correct index is need to extract the right catchment id as the catchment ids may not be in order
  this->ncats = size_v;
  for (int i=0; i<ncats; i++) { 
    this->cat_id_index[cat_id[i]-1] = i;   
  }

  
}

/*
Reads the spatial file that provides mapping between grid ids (NWM) and catchment ids and the area fraction of each grid within the catchment

*/
void smc_mapping::SoilMoistureMapping::
ReadCatchmentParamsData(std::string spatial_file)
{
  std::ifstream fp;
  fp.open(spatial_file);
  if (!fp) {
    cout<<"file "<<spatial_file<<" doesn't exist. \n";
    abort();
  }

  // vectors with _v are local vectors created to perform basic operations that we can't do with arrays directly (e.g., sorting)
  
  std::vector<double> global_deficit_v(0.0);
  std::vector<double> porosity_v(0.0);
  std::vector<double> depth_v(0.0);
  std::vector<double> szm_v(0.0);
  std::vector<string> vars;
  std::string line, cell;
  
  //read first line of strings which contains variables names.
  std::getline(fp, line);
  std::stringstream lineStream(line);
 
  while(std::getline(lineStream,cell, ',')) {
    vars.push_back(cell);
  }

  // check the order of the variables; this can be avoided
  if (vars[0].compare("global_deficit") != 0) {
    std::stringstream errMsg;
    errMsg << "data order issue: grid_id should be the first column "<<"\n";
    throw std::runtime_error(errMsg.str());
  }
  else if (vars[1].compare("porosity") != 0) {
    std::stringstream errMsg;
    errMsg << "data order issue: cat_id should be the second column "<<"\n";
    throw std::runtime_error(errMsg.str());
  }
  else if (vars[2].compare("depth") != 0) {
    std::stringstream errMsg;
    errMsg << "data order issue: area_fraction should be the third column "<<"\n";
    throw std::runtime_error(errMsg.str());
  }
  else if (vars[3].compare("szm") != 0) {
    std::stringstream errMsg;
    errMsg << "data order issue: area_fraction should be the third column "<<"\n";
    throw std::runtime_error(errMsg.str());
  }

  int count;
  
  while (std::getline(fp, line)) {
    
    std::stringstream lineStream(line);
    count = 0;
    
    while(std::getline(lineStream,cell, ',')) {
      switch (count) {
      case 0 : {
	//global_deficit_v.push_back(stod(cell));
	break;
      }
      case 1 : {
	global_deficit_v.push_back(stod(cell));
	break;
      }
      case 2 : {
	porosity_v.push_back(stod(cell));
	break;
      }
      case 3 : {
	depth_v.push_back(stod(cell));
	break;
      }
      case 4 : {
	szm_v.push_back(stod(cell));
	break;
      }
      default :
	std::cerr<<"Number of columns more than expected!";
	abort();
      }
      
      count++;
    }
    
  }
  
  int size_v = global_deficit_v.size();

  // memory allocation
  this->cat_global_deficit_V1 = new double[size_v];
  this->maxsmc_V1 = new double[size_v];
  this->depth_V1 = new double[size_v];
  this->szm_V1 = new double[size_v];

  // copyback data from vectors to 1D arrays
  for (int i=0; i<size_v; i++) {
    this->cat_global_deficit_V1[i] = global_deficit_v[i];
    this->maxsmc_V1[i] = porosity_v[i];
    this->depth_V1[i] = depth_v[i];
    this->szm_V1[i] = szm_v[i];
  }
  
}


smc_mapping::SoilMoistureMapping::
~SoilMoistureMapping()
{}

#endif
