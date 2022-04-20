#ifndef SMCP_C_INCLUDED
#define SMCP_C_INCLUDED

#include <cstring>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include "../include/soil_moisture_profile.hxx"


smc_profile::SMCProfile::
SMCProfile()
{
  this->shape[0] = 1;
  this->shape[1] = 1;
  this->shape[2] = 1;
  this->spacing[0] = 1.;
  this->spacing[1] = 1.;
  this->origin[0] = 0.;
  this->origin[1] = 0.;
  this->soil_depth =0.0;
  this->soil_storage_model= "";
  this->ncells=0;
  this->init_profile = true;
}

smc_profile::SMCProfile::
SMCProfile(std::string config_file)
{
  this->InitFromConfigFile(config_file);
  this->shape[0] = this->ncells;
  this->shape[1] = 1;
  this->shape[2] = 1;
  this->spacing[0] = 1.;
  this->spacing[1] = 1.;
  this->origin[0] = 0.;
  this->origin[1] = 0.;
  this->InitializeArrays();
  this->init_profile = true;
}

void smc_profile::SMCProfile::
InitializeArrays(void)
{
  this->soil_moisture_profile = new double[ncells];
  this->soil_moisture_layered = new double[ncells];
  this->soil_storage_m = 0.0;
  this->soil_storage_change_per_timestep_m = 0.0;
}

/*
Read and initialize values from configuration file
@input - soilZ   (1D)  :   soil discretization; array of depths from the surface [m]
@input - layersZ  (1D) : depth of each layer from the surface [m]
@input - bb  (double)  : pore size distribution [-], beta exponent on Clapp-Hornberger (1978)
@input - satpsi  (double) : saturated capillary head (saturated moisture potential) [m]
@input - ncells  (int) : number of cells of the discretized soil column
@input - nlayers (int) : numer of soil moisture layers
@input - water_table_thickness_m : thickness of the water table from the bottom of the computational domain (soil column) [m]; the config file provides the initial value. Default value is set to 0.1 m thickness
@input - soil_storage_model (string) : Conceptual or Layered soil reservoir models
@input - soil_moisture_profile_option (string) : valid only when layered soil reservoir model is chosen; option include `constant` or `linear`. The option `constant` assigns a constant value to discretized cells within each layer, `linear` option linearly interpolate values between layers and interpolated values are assigned to the soil discretization
@params - input_var_names_model (1D) : dynamically sets model inputs to be used in the bmi input_var_names
*/

void smc_profile::SMCProfile::
InitFromConfigFile(std::string config_file)
{ 
  std::ifstream fp;
  fp.open(config_file);
  
  bool is_soilZ_set = false;
  bool is_layersZ_set = false;
  bool is_smcmax_set = false;
  bool is_bb_set = false;
  bool is_satpsi_set = false;
  bool is_wt_set = false;
  bool is_soil_storage_model_set = false;
  bool is_soil_moisture_profile_option_set = false; // option for linear or piece-wise constant layered profile
  
  while (fp) {

    std::string key;
    std::getline(fp, key);
    
    int loc = key.find("=");
    std::string key_sub = key.substr(0,loc);
    
    if (key_sub == "soil_params.Z") {
      std::string tmp_key = key.substr(loc+1,key.length());
      std::vector<double> vec = ReadVectorData(tmp_key);
      
      this->soilZ = new double[vec.size()];
      
      for (unsigned int i=0; i < vec.size(); i++)
      	this->soilZ[i] = vec[i];
      
      this->ncells = vec.size();
      this->soil_depth = this->soilZ[this->ncells-1];
      is_soilZ_set = true;
      continue;
    }
    else if (key_sub == "soil_params.layersZ") {
      std::string tmp_key = key.substr(loc+1,key.length());
      std::vector<double> vec = ReadVectorData(tmp_key);
      
      this->layersZ = new double[vec.size()];
      
      for (unsigned int i=0; i < vec.size(); i++)
	this->layersZ[i] = vec[i];
      
      this->nlayers = vec.size();
      this->last_layer_depth = this->layersZ[this->nlayers-1];
      is_layersZ_set = true;
      continue;
    }
    else if (key_sub == "soil_params.smcmax") {
      this->smcmax = std::stod(key.substr(loc+1,key.length()));
      is_smcmax_set = true;
      continue;
    }
    else if (key_sub == "soil_params.b") {
      this->bb = std::stod(key.substr(loc+1,key.length()));
      assert (this->bb > 0);
      is_bb_set = true;
      continue;
    }
    else if (key_sub == "soil_params.satpsi") {
      this->satpsi = std::stod(key.substr(loc+1,key.length()));
      is_satpsi_set = true;
      continue;
    }
    else if (key_sub == "soil_params.water_table_thickness") {
      this->water_table_thickness_m = std::stod(key.substr(loc+1,key.length()));
      is_wt_set = true;
      continue;
    }
    else if (key_sub == "soil_storage_model") {
      this->soil_storage_model = key.substr(loc+1,key.length());
      is_soil_storage_model_set = true;
      continue;
    }
    else if (key_sub == "soil_moisture_profile_option") {  //Soil moisture profile option
      this->soil_moisture_profile_option = key.substr(loc+1,key.length());
      is_soil_moisture_profile_option_set = true;
      continue;
    }
  }
  fp.close();
  
  if (!is_soilZ_set) {
    std::stringstream errMsg;
    errMsg << "soilZ not set in the config file "<< config_file << "\n";
    throw std::runtime_error(errMsg.str());
  }

  if (!is_layersZ_set) {
    if (this->soil_storage_model == "layered" || this->soil_storage_model == "Layered") {
      std::stringstream errMsg;
      errMsg << "layersZ not set in the config file "<< config_file << "\n";
      throw std::runtime_error(errMsg.str());
    }
  }
  
  if (!is_smcmax_set) {
    std::stringstream errMsg;
    errMsg << "smcmax not set in the config file "<< config_file << "\n";
    throw std::runtime_error(errMsg.str());
  }
  
  if (!is_bb_set) {
    std::stringstream errMsg;
    errMsg << "bb (Clapp-Hornberger's parameter) not set in the config file "<< config_file << "\n";
    throw std::runtime_error(errMsg.str());
  }
  
  if (!is_satpsi_set) {
    std::stringstream errMsg;
    errMsg << "satpsi not set in the config file "<< config_file << "\n";
    throw std::runtime_error(errMsg.str());
  }
  
  if (!is_wt_set) {
    std::cout<<"Warning! Water table location not provided, defualt is 1.9 m deep. \n";
    //initial water table location 1.9 m deep, if not provided in the config file
    this->water_table_thickness_m = this->soil_depth - 1.9; 
  }
  
  if(is_soil_storage_model_set) {
    
    if (this->soil_storage_model == "conceptual" || this->soil_storage_model == "Conceptual") {
      input_var_names_model = new std::vector<std::string>;
      input_var_names_model->push_back("soil_storage");
      input_var_names_model->push_back("soil_storage_change");
      this->soil_moisture_profile_option_bmi = 1;

    }
    else if (this->soil_storage_model == "layered" || this->soil_storage_model == "Layered") {
      if (!is_soil_moisture_profile_option_set) {
	std::stringstream errMsg;
	errMsg << "soil moisture profile option not set in the config file "<< config_file << "\n";
	throw std::runtime_error(errMsg.str());
      }
      
      input_var_names_model = new std::vector<std::string>;
      input_var_names_model->push_back("soil_storage");
      input_var_names_model->push_back("soil_storage_change");
      input_var_names_model->push_back("soil_moisture_layered");

      this->soil_moisture_profile_option_bmi = 2;
      
    }
  }
  
  // check if the size of the input data is consistent
  assert (this->ncells > 0);
  
}

/*
returns dynamically allocated 1D vector of strings that contains correct input variable names based on the model (conceptual or layered) chosen
*/
std::vector<std::string>* smc_profile::SMCProfile::
InputVarNamesModel()
{
  return input_var_names_model;
}

/*
Reads 1D data from the config file
- used for reading soil discretization (1D)
- used for reading layers depth from the surface if model `layered` is chosen
*/
std::vector<double> smc_profile::SMCProfile::
ReadVectorData(std::string key)
{
  int pos =0;
  std::string delimiter = ",";
  std::vector<double> value(0.0);
  std::string z1 = key;
  
  if (z1.find(delimiter) == std::string::npos) {
    double v = stod(z1);
    if (v == 0.0) {
      std::stringstream errMsg;
      errMsg << "soilZ (depth of soil reservior) should be greater than zero. It it set to "<< v << " in the config file "<< "\n";
      throw std::runtime_error(errMsg.str());
    }
    
    value.push_back(v);
    
  }
  else {
    while (z1.find(delimiter) != std::string::npos) {
      pos = z1.find(delimiter);
      std::string z_v = z1.substr(0, pos);

      value.push_back(stod(z_v.c_str()));
      
      z1.erase(0, pos + delimiter.length());
      if (z1.find(delimiter) == std::string::npos)
	value.push_back(stod(z1));
    }
  }
  
  return value;
}

/*
- Computes 1D soil moisture profile for conceptual reservoir using Newton-Raphson iterative method
- local_variables:
  @param lam  [-] : 1/bb (bb: pore size distribution)
  @param satpsi_cm [cm] : saturated moisture potential
  @param depth  [cm] : depth of the soil column
  @param zb  [cm] : bottom of the computational domain
  @param z0  [cm] : bottom of the fictitious domain (to track fictitious water table location)
  @param z1  [cm] : depth of the water table at the previous timestep (known)
  @param z2  [cm] : depth of the water table at the current timestep (unknown, will use the Newton-Raphson method to find it)
  @param soil_storage_max [cm] : maximum soil storage
  @param soil_storage_previous_timestepcm [cm] : soil storage at the previous timestep computed locally from the watertable location at the previous timestep
  @param soil_storage_current_timestepcm [cm] : soil storage at the current timestep
  @param tol [cm] : Error tolerance for finding the new root (i.e., water_table_thickness)
  @param soil_moisture_profile [-] : OUTPUT (soil moisture content vertical profile [-])
*/
void smc_profile::SMCProfile::
SoilMoistureProfileFromConceptualReservoir()
{
  // converting variables to cm for numerical reasons only
  double satpsi_cm = this->satpsi * 100.;
  double depth = this->soil_depth * 100.;
  double zb=0; // bottom of the computational domain
  double z0=0; // bottom of the fictitious domain (to track fictitious water table location)
  
  double zi = 0.01; // initial guess for the water table location, use Newton-Raphson to find new zi
  double soil_storage_max = depth * this->smcmax;
  
  double soil_storage_change_per_timestep_cm = soil_storage_change_per_timestep_m * 100.0;
  double soil_storage_current_timestep_cm = 100.0 * this->soil_storage_m;  /* storage at the current timestep */

  double lam = 1.0/this->bb;
  double beta = 1.0 - lam;
  double alpha = pow(satpsi_cm,lam)/beta; // a constant term obtained in the integration of the soil moisture function
  double tol = 0.000001;

  int count = 0;

  // compute a new profile only if sufficient amount of water is added at this timestep. 1.0E-4 corresponds to 0.001 mm of water
  
  if (soil_storage_change_per_timestep_cm > 0.0001 || this->init_profile) {
    
    // turn of the flag for times t > 0 
    this->init_profile = false;
    
    // check if the storage is greater than the maximum soil storage. if yes, set it to the maximum storage
    if(soil_storage_current_timestep_cm >= soil_storage_max) {
      for(int j=0;j<this->ncells;j++)
	this->soil_moisture_profile[j] = this->smcmax;
      return;
    }
    
    double diff=1000.0; // guess for the initial differnce between the roots
  
    double f,  zi_new, df_dzi;
    
    do {
      count++;
      
      if(count>10000) {
	throw std::runtime_error("No convergence loop count: after 10000 iterations!");
      }
      
      // function representing the total amount of soil moisture. 2nd term is the integral of the Clap-Hornberger function (area under the soil moisture curve)
      //f = this->smcmax * (z2-z1) + alpha * this->smcmax * (pow(depth-z2,beta)- pow(depth-z1,beta)) -  soil_storage_change_per_timestep_cm;

      double fis = this->smcmax * (zi - z0) + this->smcmax * satpsi_cm + alpha * this->smcmax * ( pow((depth-zi),beta) - pow(satpsi_cm,beta) );
      double fib = this->smcmax * (zi - z0) + this->smcmax * satpsi_cm + alpha * this->smcmax * ( pow(abs(zb-zi),beta) - pow(satpsi_cm,beta) );
	

      fib = zi >= 0.0 ? 0.0 : fib;

      f = fis - fib - soil_storage_current_timestep_cm;

      // derivative of f w.r.t zi
      double dfis = this->smcmax - alpha * this->smcmax * beta * pow((depth-zi),beta-1.0);
      double dfib = this->smcmax - alpha * this->smcmax * beta * pow(std::fabs(zb-zi),beta-1.0);

      dfib = zi >= 0.0 ? 0.0 : dfib;
      df_dzi = dfis - dfib;
      
      // Newton-Raphson method
      zi_new = zi - f / std::max(df_dzi,1.e-6); // to avoid division by zero
      
      diff=zi_new-zi; // difference betweent the previous and new root
      
      zi=zi_new;     // update the previous root

      z0 = zi >= 0.0 ? zb : zi - satpsi_cm;
      
	
    } while (std::fabs(diff) < tol);

    
    this->water_table_thickness_m = zi/100.;
    
    
    /*******************************************************************/
    /* get a high resolution moisture profile that will be mapped on the desired soil discretization */
    
    int z_hres = 1000;
    double *smct_temp = new double[z_hres];
    double *z_temp = new double[z_hres];
    
    // we have the new water table location now, so let's compute the soil moisture curve now
    double z = zi + satpsi_cm;
    double dz_v = (depth - zi - satpsi_cm)/z_hres; // vertical spacing
    double z_head = satpsi_cm; // input variable to the soil moisture function
    
    for (int i=0;i<z_hres;i++) {
      z_head += dz_v;
      smct_temp[i] = this->smcmax * std::pow((satpsi_cm/z_head),lam) ;
      
      z+=dz_v;
      z_temp[i] = z;
    }
    
    // map the high resolution soil moisture curve to the soil discretization depth that is provided in the config file
    for (int i=0; i<this->ncells; i++) {
      for (int j=0; j<z_hres; j++) {
	if ( depth - soilZ[ncells-1-i]*100 <= zi + satpsi_cm)
	  this->soil_moisture_profile[i] = smcmax;
	else if (z_temp[j]  >= (depth - this->soilZ[i]*100) ) {
	  this->soil_moisture_profile[i] = smct_temp[j];
	  break;
	}
      }
    }
    
  }
  
}

/*
- Computes 1D soil moisture profile for layered/calculated reservoir. That is, take layered reservior soil moisture and distribute it vertically to the desired soil discretization
- Two strategies are implemented
  - Constant strategy: A simple technique to map layered values to grids in the soil discretization. That is, all grid cells in the discretization has a constant value within a layer. Note cells at the interface take average value of the layers 
  - Linear strategy: A linear interpolation tehcnique is used to map layered values to grids in the soil discretization. That is, grid cells in the discretization take linearly interpolated value between consecutive layers.
 - Note: the water table location is the thickness of the water table plus saturated capillary head (satpsi)
- local_variables:
  @param lam  [-]               : 1/bb (bb: pore size distribution)
  @param soil_depth  [cm]       : depth of the soil column
  @param last_layer_depth  [cm] : depth of the last layer from the surface
*/

void smc_profile::SMCProfile::
SoilMoistureProfileFromLayeredReservoir()
{
  double lam=1.0/this->bb; // pore distribution index
  double water_table_thickness = soil_depth - this->satpsi - this->water_table_thickness_m; // water table depth at the current timestep
  std::vector<double> z_layers_n(1,0.0);
  
  for (int i=0; i <this->nlayers; i++)
    z_layers_n.push_back(layersZ[i]);

  std::vector<double> smc_column;
  int c = 0;
  double delta = 0.0;

  // piece-wise constant (vertically)
  if (this->soil_moisture_profile_option == "constant" || this->soil_moisture_profile_option == "Constant") {
    bool layers_flag=true;
    
    // loop over all the cells in the discretized column
    for (int i=0; i < ncells; i++) {
      
      if (soilZ[i] < layersZ[c]) {  // cell completely lie within a layer
	
	this->soil_moisture_profile[i] = soil_moisture_layered[c];
	
      }
      else if (soilZ[i] < this->last_layer_depth) { // cell at the interface of layers, so take the mean
	
	this->soil_moisture_profile[i] = 0.5*(soil_moisture_layered[c] + soil_moisture_layered[c+1]);
	c++;
	
      }
      else { // extend the profile below the last layer depth
	
	double zz = this->soil_depth - soilZ[i];
	double theta = delta + std::pow((this->satpsi/zz),lam)*this->smcmax;
	
	if (layers_flag) { // Ensure continuity of soil moisture at the interface of the depth of the last layer and profile below that depth of the last layer
	  delta = 0.0 - theta;
	  theta = theta + delta;
	  layers_flag=false;
	}
	
	double theta1 = soil_moisture_layered[this->nlayers-1] +  theta;
	double theta2 = std::min(this->smcmax, theta1);

	this->soil_moisture_profile[i] = soilZ[i] <  water_table_thickness ? theta2 : this->smcmax;
	/*
	if (soilZ[i] < soil_depth - this->satpsi - this->water_table_depth_m)
	  this->soil_moisture_profile[i] = theta2;
	else
	this->soil_moisture_profile[i] = this->smcmax;*/
      }
      
    }
    
  }
  else if (this->soil_moisture_profile_option == "linear" || this->soil_moisture_profile_option == "Linear") {
    
    bool layers_flag=true;
    double t_v=0.0;
   
    this->soil_moisture_profile[0] = soil_moisture_layered[0]; // first cell gets the top-layer soil moisture

    for (int i=1; i < ncells; i++) {

      // linear interpolation between consecutive layers
      if (soilZ[i] <= layersZ[c] && c < this->nlayers-1) {
	t_v = LinearInterpolation(z_layers_n[c], z_layers_n[c+1], soil_moisture_layered[c], soil_moisture_layered[c+1], soilZ[i]);
	this->soil_moisture_profile[i] = t_v;
	
	if (soilZ[i+1] > layersZ[c]) // the parameter c keeps track of the layer
	  c++;
      }
      else {
	// extend the profile below the depth of the last layer
	double zz = this->soil_depth - soilZ[i];
	double theta = delta + std::pow((this->satpsi/zz),lam)*this->smcmax;
	
	if (layers_flag) { // Ensure continuity of soil moisture at the interface of the depth of the last layer and profile below that depth of the last layer
	  delta = 0.0 - theta;
	  theta = theta + delta;
	  layers_flag = false;
	}
	
	double theta1 = soil_moisture_layered[this->nlayers-1] +  theta;
	double theta2 = std::min(this->smcmax, theta1);

	this->soil_moisture_profile[i] = soilZ[i] <  water_table_thickness ? theta2 : this->smcmax;
	
      }
    }

  }
  else {
    std::stringstream errMsg;
    errMsg << "Soil moisture profile "<< this->soil_storage_model << " works with options \'constant\' and \'linear\'. Provide at least one option in the config file "<<"\n";
    throw std::runtime_error(errMsg.str());
  }

  // for vis comparison with the python version
  for (int j=0; j< ncells; j++)
    std::cout<<soilZ[j] <<" "<<this->soil_moisture_profile[j]<<"\n";
  
  // find and update water table location
  for (int j=0; j< ncells; j++) {
    if (this->soil_moisture_profile[j] == this->smcmax) {
      this->water_table_thickness_m = soil_depth - soilZ[j] - this->satpsi; // check with Fred??
      break;
    }
  }
  
}

double smc_profile::SMCProfile::
LinearInterpolation(double z1, double z2, double t1, double t2, double z)
{

  double m = (t2 - t1) / (z2 - z1);
  return t1 + m * (z - z1);

}

smc_profile::SMCProfile::
~SMCProfile()
{}

#endif
