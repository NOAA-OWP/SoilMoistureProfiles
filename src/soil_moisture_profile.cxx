#ifndef SMP_C_INCLUDED
#define SMP_C_INCLUDED

#include <cstring>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include "../include/soil_moisture_profile.hxx"

enum {Conceptual=1, Layered=2};

soil_moisture_profile::SoilMoistureProfile::
SoilMoistureProfile()
{
  this->shape[0] = 1;
  this->shape[1] = 1;
  this->shape[2] = 1;
  this->spacing[0] = 1.;
  this->spacing[1] = 1.;
  this->origin[0] = 0.;
  this->origin[1] = 0.;
  this->soil_depth =0.0;
  this->soil_storage_model= 1;
  this->ncells=0;
  this->init_profile = true;
}

soil_moisture_profile::SoilMoistureProfile::
SoilMoistureProfile(std::string config_file)
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

void soil_moisture_profile::SoilMoistureProfile::
InitializeArrays(void)
{
  this->soil_moisture_profile = new double[ncells];
  this->soil_moisture_layered = new double[ncells];
  this->soil_storage = 0.0;
  this->soil_storage_change_per_timestep = 0.0;
}

/*
Read and initialize values from configuration file
@input - soil_z   (1D)  :   soil discretization; array of depths from the surface [m]
@input - layers_z  (1D) : depth of each layer from the surface [m]
@input - bb  (double)  : pore size distribution [-], beta exponent on Clapp-Hornberger (1978)
@input - satpsi  (double) : saturated capillary head (saturated moisture potential) [m]
@input - ncells  (int) : number of cells of the discretized soil column
@input - nlayers (int) : numer of soil moisture layers
@input - soil_storage_model (string) : Conceptual or Layered soil reservoir models
@input - soil_moisture_profile_option (string) : valid only when layered soil reservoir model is chosen; option include `constant` or `linear`. The option `constant` assigns a constant value to discretized cells within each layer, `linear` option linearly interpolate values between layers and interpolated values are assigned to the soil discretization
@params - input_var_names_model (1D) : dynamically sets model inputs to be used in the bmi input_var_names
*/

void soil_moisture_profile::SoilMoistureProfile::
InitFromConfigFile(std::string config_file)
{ 
  std::ifstream fp;
  fp.open(config_file);
  
  bool is_soil_z_set = false;
  bool is_layers_z_set = false;
  bool is_smcmax_set = false;
  bool is_bb_set = false;
  bool is_satpsi_set = false;
  bool is_soil_storage_model_set = false;
  bool is_soil_moisture_layered_option_set = false; // option for linear or piece-wise constant layered profile
  bool is_soil_storage_model_depth_set = false;
  
  while (fp) {

    std::string line;
    std::string param_key, param_value, param_unit;
    
    std::getline(fp, line);
   
    int loc_eq = line.find("=") + 1;
    int loc_u = line.find("[");
    param_key = line.substr(0,line.find("="));

    bool is_unit = line.find("[") != std::string::npos;

    if (is_unit)
      param_unit = line.substr(loc_u,line.find("]")+1);
    else
      param_unit = "";

    param_value = line.substr(loc_eq,loc_u - loc_eq);
    
    if (param_key == "soil_z") {
      std::vector<double> vec = ReadVectorData(param_value);
      
      this->soil_z = new double[vec.size()];
      
      for (unsigned int i=0; i < vec.size(); i++)
      	this->soil_z[i] = vec[i];
      
      this->ncells = vec.size();
      this->soil_depth = this->soil_z[this->ncells-1];
      is_soil_z_set = true;
      continue;
    }
    else if (param_key == "soil_layers_z") {
      std::vector<double> vec = ReadVectorData(param_value);
      this->layers_z = new double[vec.size()];
      
      for (unsigned int i=0; i < vec.size(); i++)
	this->layers_z[i] = vec[i];
      
      this->nlayers = vec.size();
      this->last_layer_depth = this->layers_z[this->nlayers-1];
      is_layers_z_set = true;
      continue;
    }
    else if (param_key == "soil_params.smcmax") {
      this->smcmax = std::stod(param_value);
      is_smcmax_set = true;
      continue;
    }
    else if (param_key == "soil_params.b") {
      this->bb = std::stod(param_value);
      assert (this->bb > 0);
      is_bb_set = true;
      continue;
    }
    else if (param_key == "soil_params.satpsi") {
      this->satpsi = std::stod(param_value);
      is_satpsi_set = true;
      continue;
    }
    else if (param_key == "soil_storage_model") {
      if ( param_value == "Conceptual" || param_value == "conceptual")
	this->soil_storage_model = Conceptual;
      else if (param_value == "layered" || param_value == "Layered") 
	this->soil_storage_model = Layered;

      is_soil_storage_model_set = true;
      continue;
    }
    else if (param_key == "soil_moisture_layered_option") {  //Soil moisture profile option
      this->soil_moisture_layered_option = stod(param_value);
      is_soil_moisture_layered_option_set = true;
      continue;
    }
    else if (param_key == "soil_storage_depth") {
      this->soil_storage_model_depth = std::stod(param_value);
      assert (this->soil_storage_model_depth > 0);
      is_soil_storage_model_depth_set = true;
      continue;
    }
  }
  
  fp.close();
  
  if (!is_soil_z_set) {
    std::stringstream errMsg;
    errMsg << "soil_z not set in the config file "<< config_file << "\n";
    throw std::runtime_error(errMsg.str());
  }

  if (!is_layers_z_set) {
    if (this->soil_storage_model == Layered) {
      std::stringstream errMsg;
      errMsg << "layers_z not set in the config file "<< config_file << "\n";
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

  if (!is_soil_storage_model_depth_set) {
    std::stringstream errMsg;
    errMsg << "soil_storage_model_depth not set in the config file "<< config_file << "\n";
    throw std::runtime_error(errMsg.str());
  }

    
  if(is_soil_storage_model_set) {
    
    if (this->soil_storage_model == Conceptual) {
      input_var_names_model = new std::vector<std::string>;
      input_var_names_model->push_back("soil_storage");
      input_var_names_model->push_back("soil_storage_change");
    }
    else if (this->soil_storage_model == Layered) {
      if (!is_soil_moisture_layered_option_set) {
	std::stringstream errMsg;
	errMsg << "soil moisture profile option not set in the config file "<< config_file << "\n";
	throw std::runtime_error(errMsg.str());
      }
      
      input_var_names_model = new std::vector<std::string>;
      input_var_names_model->push_back("soil_storage");
      input_var_names_model->push_back("soil_storage_change");
      input_var_names_model->push_back("soil_moisture_layered");
    }
  }
  
  // check if the size of the input data is consistent
  assert (this->ncells > 0);
  
}

/*
returns dynamically allocated 1D vector of strings that contains correct input variable names based on the model (conceptual or layered) chosen
*/
std::vector<std::string>* soil_moisture_profile::SoilMoistureProfile::
InputVarNamesModel()
{
  return input_var_names_model;
}

/*
Reads 1D data from the config file
- used for reading soil discretization (1D)
- used for reading layers depth from the surface if model `layered` is chosen
*/
std::vector<double> soil_moisture_profile::SoilMoistureProfile::
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
      errMsg << "soil_z (depth of soil reservior) should be greater than zero. It it set to "<< v << " in the config file "<< "\n";
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


void soil_moisture_profile::SoilMoistureProfile::
SoilMoistureProfileUpdate()
{
  if (this->soil_storage_model == Conceptual) {
    this->SoilMoistureProfileFromConceptualReservoir();
  }
  else if (this->soil_storage_model == Layered) {
    this->SoilMoistureProfileFromLayeredReservoir();
  }
  else {
    std::stringstream errMsg;
    errMsg << "Soil moisture profile OPTION provided in the config file is " << this->soil_storage_model<< ", which should be either \'concepttual\' or \'layered\' " <<"\n";
    throw std::runtime_error(errMsg.str());
  }

}
/*
  Computes 1D soil moisture profile for conceptual reservoir using Newton-Raphson iterative method
  For detailed decription of the model implemented here, please see README.md on the github repo
  local_variables:
  @param lam  [-] : 1/bb (bb: pore size distribution)
  @param satpsi_cm [cm] : saturated moisture potential
  @param soil_depth_cm  [cm] : depth of the soil column (note this depth can be different than the depth of the soil_storage_model (e.g., CFE)
  @param soil_storage_model_depth [m] : depth of the soil reservoir model (e.g., CFE)
  @param zb  [cm] : bottom of the computational domain
  @param z0  [cm] : bottom of the fictitious domain (to track fictitious water table location)
  @param z1  [cm] : depth of the water table at the previous timestep (known)
  @param z2  [cm] : depth of the water table at the current timestep (unknown, will use the Newton-Raphson method to find it)
  @param soil_storage_max [cm] : maximum soil storage
  @param soil_storage_previous_timestepcm [cm] : soil storage at the previous timestep computed locally from the watertable location at the previous timestep
  @param soil_storage_current_timestepcm [cm] : soil storage at the current timestep
  @param tol [cm] : Error tolerance for finding the new root (i.e., water_table_thickness)
  @param soil_moisture_profile [-] : OUTPUT (soil moisture content vertical profile [-])
  ** NOTE: the module needs to be fixed if the CFE and SFT soil depths are different
*/
void soil_moisture_profile::SoilMoistureProfile::
SoilMoistureProfileFromConceptualReservoir()
{
  // converting variables to cm for numerical reasons only
  double satpsi_cm = this->satpsi * 100.;
  //double soil_depth_cm = this->soil_depth * 100.; // soil profile depth
  double model_depth = this->soil_storage_model_depth * 100.;
  double zb=0; // bottom of the computational domain
  double z0=0; // bottom of the fictitious domain (to track fictitious water table location)
  
  double zi = 0.01; // initial guess for the water table location, use Newton-Raphson to find new zi
  double soil_storage_max = model_depth * this->smcmax;
  
  double soil_storage_change_per_timestep_cm = soil_storage_change_per_timestep * 100.0;
  double soil_storage_current_timestep_cm = 100.0 * this->soil_storage;  /* storage at the current timestep */

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
  
    double f, zi_new, df_dzi;
    
    do {
      count++;
      
      if(count>10000) {
	throw std::runtime_error("No convergence loop count: after 10000 iterations!");
      }
      
      // function representing the total amount of soil moisture. 2nd term is the integral of the Clap-Hornberger function (area under the soil moisture curve)

      // fis is integrates the function from z0 to the surface: 1st part: saturated soil between zi and z0; 2nd part: capillary fringe; 3rd: enclosed area between satpis and the surface
      // fib is non-zero for zi <0, zero otherwise. this is the volume of water that needs to be subtracted from "fis" to get the water storage in the computational domain (say top 2 m if soil column has depth 2m)
      double fis = this->smcmax * (zi - z0) + this->smcmax * satpsi_cm + alpha * this->smcmax * ( pow((model_depth-zi),beta) - pow(satpsi_cm,beta) );
      double fib = this->smcmax * (zi - z0) + this->smcmax * satpsi_cm + alpha * this->smcmax * ( pow(abs(zb-zi),beta) - pow(satpsi_cm,beta) );
	

      fib = zi >= 0.0 ? 0.0 : fib;

      f = fis - fib - soil_storage_current_timestep_cm;

      // derivative of f w.r.t zi
      double dfis = this->smcmax - alpha * this->smcmax * beta * pow((model_depth-zi),beta-1.0);
      double dfib = this->smcmax - alpha * this->smcmax * beta * pow(fabs(zb-zi),beta-1.0);

      dfib = zi >= 0.0 ? 0.0 : dfib;
      df_dzi = dfis - dfib;
      
      // Newton-Raphson method
      zi_new = zi - f / max(df_dzi,1.e-6); // to avoid division by zero
      
      diff=zi_new-zi; // difference betweent the previous and new root
      
      zi=zi_new;     // update the previous root

      z0 = zi >= 0.0 ? zb : zi - satpsi_cm;
      
      //std::cout<<"water table: "<<count<<" "<<zi <<" "<<fis<<" "<<fib<<" "<<f<<" "<<dfis<<" "<<dfib<<" "<<df_dzi<<": "<<diff<<" "<<std::fabs(diff)<<" "<<tol<<"\n";	
    } while (fabs(diff) > tol);

    // water table thickness can be negative and that would be depth of the water table below the depth of the computational domain; probably a better name would be water_table_location
    this->water_table_thickness = zi/100.;
    
    /*******************************************************************/
    /* get a high resolution moisture profile that will be mapped on the desired soil discretization */
    
    int z_hres = 1000;
    double *smct_temp = new double[z_hres];
    double *z_temp = new double[z_hres];
    
    // we have the new water table location now, so let's compute the soil moisture curve now
    double z = zi + satpsi_cm;
    double dz_v = (model_depth - zi - satpsi_cm)/z_hres; // vertical spacing
    double z_head = satpsi_cm; // input variable to the soil moisture function
    
    for (int i=0;i<z_hres;i++) {
      z_head += dz_v;
      smct_temp[i] = this->smcmax * pow((satpsi_cm/z_head),lam) ;
      
      z+=dz_v;
      z_temp[i] = z;
    }
    
    // map the high resolution soil moisture curve to the soil discretization depth that is provided in the config file
    for (int i=0; i<this->ncells; i++) {
      for (int j=0; j<z_hres; j++) {
	if ( (model_depth - soil_z[i]*100) <= (zi + satpsi_cm) )
	  this->soil_moisture_profile[i] = smcmax;
	else if (z_temp[j]  >= (model_depth - this->soil_z[i]*100) ) {
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

void soil_moisture_profile::SoilMoistureProfile::
SoilMoistureProfileFromLayeredReservoir()
{
  double lam=1.0/this->bb; // pore distribution index

  //double water_table_thickness = this->soil_depth - this->satpsi - this->water_table_thickness_m; // water table depth at the current timestep -- this is not needed??
  
  std::vector<double> z_layers_n(1,0.0);
  
  for (int i=0; i <this->nlayers; i++)
    z_layers_n.push_back(layers_z[i]);

  std::vector<double> smc_column;
  int c = 0;
  double delta = 0.0;

  // piece-wise constant (vertically)
  if (this->soil_moisture_layered_option == "constant" || this->soil_moisture_layered_option == "Constant") {
    bool layers_flag=true;
    
    // loop over all the cells in the discretized column
    for (int i=0; i < ncells; i++) {
      
      if (soil_z[i] < layers_z[c]) {  // cell completely lie within a layer
	
	this->soil_moisture_profile[i] = soil_moisture_layered[c];
	
      }
      else if (soil_z[i] < this->last_layer_depth) { // cell at the interface of layers, so take the mean
	
	this->soil_moisture_profile[i] = 0.5*(soil_moisture_layered[c] + soil_moisture_layered[c+1]);
	c++;
	
      }
      else { // extend the profile below the last layer depth
	
	double zz = this->soil_depth - soil_z[i];
	double theta = delta + std::pow((this->satpsi/zz),lam)*this->smcmax;
	
	if (layers_flag) { // Ensure continuity of soil moisture at the interface of the depth of the last layer and profile below that depth of the last layer
	  delta = 0.0 - theta;
	  theta = theta + delta;
	  layers_flag=false;
	}
	
	double theta1 = soil_moisture_layered[this->nlayers-1] +  theta;
	/*
	double theta2 = std::min(this->smcmax, theta1);
	this->soil_moisture_profile[i] = soil_z[i] <  water_table_thickness ? theta2 : this->smcmax;
	*/
	
	this->soil_moisture_profile[i] = std::min(this->smcmax, theta1);


	/*
	if (soil_z[i] < soil_depth - this->satpsi - this->water_table_depth_m)
	  this->soil_moisture_profile[i] = theta2;
	else
	this->soil_moisture_profile[i] = this->smcmax;*/
      }
      
    }
    
  }
  else if (this->soil_moisture_layered_option == "linear" || this->soil_moisture_layered_option == "Linear") {
    
    bool layers_flag=true;
    double t_v=0.0;
   
    this->soil_moisture_profile[0] = soil_moisture_layered[0]; // first cell gets the top-layer soil moisture

    for (int i=1; i < ncells; i++) {

      // linear interpolation between consecutive layers
      if (soil_z[i] <= layers_z[c] && c < this->nlayers-1) {
	t_v = LinearInterpolation(z_layers_n[c], z_layers_n[c+1], soil_moisture_layered[c], soil_moisture_layered[c+1], soil_z[i]);
	this->soil_moisture_profile[i] = t_v;
	
	if (soil_z[i+1] > layers_z[c]) // the parameter c keeps track of the layer
	  c++;
      }
      else {
	// extend the profile below the depth of the last layer
	double zz = this->soil_depth - soil_z[i];
	double theta = delta + std::pow((this->satpsi/zz),lam)*this->smcmax;
	
	if (layers_flag) { // Ensure continuity of soil moisture at the interface of the depth of the last layer and profile below that depth of the last layer
	  delta = 0.0 - theta;
	  theta = theta + delta;
	  layers_flag = false;
	}
	
	double theta1 = soil_moisture_layered[this->nlayers-1] +  theta;
	/*
	double theta2 = std::min(this->smcmax, theta1);
	this->soil_moisture_profile[i] = soil_z[i] <  water_table_thickness ? theta2 : this->smcmax;
	*/

	this->soil_moisture_profile[i] = std::min(this->smcmax, theta1);
	
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
    std::cout<<soil_z[j] <<" "<<this->soil_moisture_profile[j]<<"\n";
  /*
  // find and update water table location
  for (int j=0; j< ncells; j++) {
    if (this->soil_moisture_profile[j] == this->smcmax) {
      this->water_table_thickness = soil_depth - soil_z[j] - this->satpsi; // check with Fred?? do we actually care about the location, all we need is the mass of water to be consistent
      break;
    }
  }
  */
}

double soil_moisture_profile::SoilMoistureProfile::
LinearInterpolation(double z1, double z2, double t1, double t2, double z)
{

  double m = (t2 - t1) / (z2 - z1);
  return t1 + m * (z - z1);

}

soil_moisture_profile::SoilMoistureProfile::
~SoilMoistureProfile()
{}

#endif
