#ifndef SMP_C_INCLUDED
#define SMP_C_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cmath> // fmin and fmax (C)
//#include <algorithm>  // for std::min and std::max (C++)

#include "../include/soil_moisture_profile.hxx"

enum {Conceptual=1, Layered=2, Topmodel=3};
enum {Constant=1, Linear=2};
enum {Flux_based=1, Deficit_based=2};

// soil_moisture_profile is the namespacing


void soil_moisture_profile::
SoilMoistureProfile(string config_file, struct soil_profile_parameters* parameters)
{
  
  InitFromConfigFile(config_file, parameters);

  parameters->shape[0] = parameters->ncells; // this is used for the size of soil_moisture_profile (bmi output)
  
  if (parameters->soil_storage_model == Conceptual || parameters->soil_storage_model == Topmodel)
    parameters->shape[1] = 1;
  else if (parameters->soil_storage_model == Layered)
    parameters->shape[1] = parameters->num_wetting_fronts;

  parameters->soil_moisture_profile = new double[parameters->ncells];

  parameters->soil_moisture_wetting_fronts = new double[parameters->shape[1]]();
  parameters->soil_depth_wetting_fronts    = new double[parameters->shape[1]]();

  // For water_table_based_method
  parameters->cat_area     = 1.0;        // catchment area used in the topmodel (normalized)
  parameters->shape[2]     = 1;
  parameters->spacing[0]   = 1.;
  parameters->spacing[1]   = 1.;
  parameters->origin[0]    = 0.;
  parameters->origin[1]    = 0.;
  parameters->soil_storage = 0.0;
  parameters->init_profile = true;
  parameters->soil_depth_NWM = 2.0;
  parameters->num_wetting_fronts = parameters->shape[1];
  parameters->soil_storage_change_per_timestep = 0.0;
}

/*
  Read and initialize values from configuration file
  @input soil_z   (1D)        :   soil discretization; array of depths from the surface [m]
  @input layers_z  (1D)       : depth of each layer from the surface [m]
  @input b  (double)          : pore size distribution [-], beta exponent on Clapp-Hornberger (1978)
  @input satpsi  (double)     : saturated capillary head (saturated moisture potential) [m]
  @input ncells  (int)        : number of cells of the discretized soil column
  @input ncells_layered       : number of layers (wetting fronts) for the layered model
  @input max_ncells_layered   : maximum number of layers (wetting fronts)
  @input soil_storage_model (string)           : Conceptual or Layered soil reservoir models
  @input soil_moisture_profile_option (string) : valid only when layered soil reservoir model is chosen;
                                                 option include `constant` or `linear`. The option `constant`
						 assigns a constant value to discretized cells within each layer,
						 `linear` option linearly interpolate values between layers and
						  interpolated values are assigned to the soil discretization
  @params input_var_names_model (1D)      : dynamically sets model inputs to be used in the bmi input_var_names
  @input water_table_depth (double)       : water table depth, default is 6 m for layered model, conceptual model
                                            computes its own
  @param soil_moisture_fraction     [-]   : fraction of soil moisture in the top 40 cm or user-defined
                                            soil_moisture_fraction_depth (water in the topsoil over total water)
  @param soil_moisture_fraction_depth [m] : user specified depth for the fraction of soil moisture (default is 40 cm)
  @param water_table_based_method    (string) : Flux_based (for iteration-based water table)
                                                Deficit_based (for deficit-based water table) and only needs to be provided
						if soil_storage_model is Topmodel
*/


void soil_moisture_profile::
InitFromConfigFile(string config_file, struct soil_profile_parameters* parameters)
{ 
  ifstream fp;
  fp.open(config_file);
  
  bool is_soil_z_set                    = false;
  bool is_soil_depth_layers_set         = false;
  bool is_smcmax_set                    = false;
  bool is_b_set                        = false;
  bool is_satpsi_set                    = false;
  bool is_soil_storage_model_set        = false;
  bool is_soil_storage_model_depth_set  = false;
  bool is_water_table_depth_set         = false;
  bool is_water_table_based_method_set  = false;
  bool is_soil_moisture_fraction_depth_set        = false;
  bool is_soil_moisture_profile_option_set = false; // option for linear or piece-wise constant layered profile
  
  while (fp) {

    string line;
    string param_key, param_value, param_unit;
    
    getline(fp, line);
   
    int loc_eq = line.find("=") + 1;
    int loc_u = line.find("[");
    param_key = line.substr(0,line.find("="));

    bool is_unit = line.find("[") != string::npos;

    if (is_unit)
      param_unit = line.substr(loc_u,line.find("]")+1);
    else
      param_unit = "";

    param_value = line.substr(loc_eq,loc_u - loc_eq);
    
    if (param_key == "soil_z") {
      vector<double> vec = ReadVectorData(param_key, param_value);
      
      parameters->soil_z = new double[vec.size()];
      
      for (unsigned int i=0; i < vec.size(); i++)
      	parameters->soil_z[i] = vec[i];
      
      parameters->ncells = vec.size();
      parameters->soil_depth = parameters->soil_z[parameters->ncells-1];
      is_soil_z_set = true;
      continue;
    }
    else if (param_key == "soil_depth_layers") {
      if (param_value == "bmi" || param_value == "BMI") {
	parameters->soil_depth_layers_bmi = true;
      }
      else {
	vector<double> vec = ReadVectorData(param_key,param_value);
	parameters->soil_depth_layers = new double[vec.size()];
	
	for (unsigned int i=0; i < vec.size(); i++)
	  parameters->soil_depth_layers[i] = vec[i];
	
	parameters->num_layers         = vec.size();
	parameters->num_wetting_fronts = vec.size();
	parameters->last_layer_depth   = parameters->soil_depth_layers[parameters->num_layers - 1];
	is_soil_depth_layers_set       = true;
	parameters->soil_depth_layers_bmi = false;
      }
      continue;
    }
	// NOTE: `soil_params.smcmax` may be deprecated in the future in favor of `smcmax`
    else if (param_key == "smcmax" || param_key == "soil_params.smcmax") {
      if (param_value == "bmi" || param_value == "BMI") {
	parameters->smcmax_bmi = true;
      }
      else {
	vector<double> vec = ReadVectorData(param_key, param_value);
	parameters->smcmax = new double[vec.size()];
	
	for (unsigned int i=0; i < vec.size(); i++) {
	  assert (vec[i] > 0);
	  parameters->smcmax[i] = vec[i];
	}

	parameters->smcmax_bmi = false;
	parameters->num_layers = vec.size();
	assert (parameters->num_layers > 0);
	is_smcmax_set = true;
      }
      
      continue;
    }
	// NOTE: `soil_params.b` may be deprecated in the future in favor of `b`
    else if (param_key == "b" || param_key == "soil_params.b") {
      parameters->b = stod(param_value);
      assert (parameters->b > 0);
      is_b_set = true;
      continue;
    }
	// NOTE: `soil_params.satpsi` may be deprecated in the future in favor of `satpsi`
    else if (param_key == "satpsi" || param_key == "soil_params.satpsi") {
      parameters->satpsi = stod(param_value);
      assert (parameters->satpsi > 0.0);
      is_satpsi_set = true;
      continue;
    }
    else if (param_key == "soil_storage_model") {
      if ( param_value == "Conceptual" || param_value == "conceptual")
	parameters->soil_storage_model = Conceptual;
      else if (param_value == "layered" || param_value == "Layered") 
	parameters->soil_storage_model = Layered;
      else if (param_value == "topmodel" || param_value == "TopModel" || param_value == "TOPMODEL") 
	parameters->soil_storage_model = Topmodel;

      is_soil_storage_model_set = true;
      continue;
    }
    else if (param_key == "soil_moisture_profile_option") {  //Soil moisture profile option
      if (param_value == "Constant" || param_value == "constant")
	parameters->soil_moisture_profile_option = Constant;
      else if (param_value == "Linear" || param_value == "linear")
	parameters->soil_moisture_profile_option = Linear;
      is_soil_moisture_profile_option_set = true;
      continue;
    }
    else if (param_key == "soil_storage_depth") {
      parameters->soil_storage_model_depth = stod(param_value);
      assert (parameters->soil_storage_model_depth > 0);
      is_soil_storage_model_depth_set = true;
      continue;
    }
    else if (param_key == "water_table_depth") {
      parameters->water_table_depth = stod(param_value);
      is_water_table_depth_set = true;
      continue;
    }
    else if (param_key == "soil_moisture_fraction_depth") {
      parameters->soil_moisture_fraction_depth = stod(param_value);
      is_soil_moisture_fraction_depth_set = true;
      continue;
    }
    else if (param_key == "water_table_based_method") {
      if (param_value == "deficit_based")
	parameters->water_table_based_method = Deficit_based;
      else if (param_value == "flux_based")
	parameters->water_table_based_method = Flux_based;
      is_water_table_based_method_set = true;
      continue;
    }
    else if (param_key == "verbosity") {
      
      if (param_value == "high" || param_value == "low")
	parameters->verbosity = param_value;
      else
	parameters->verbosity = "none";
	  
      continue;
    }
    
  }

  
  fp.close();
  
  if (!is_soil_z_set) {
    stringstream errMsg;
    errMsg << "soil_z not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }


  if (!is_soil_depth_layers_set && !parameters->soil_depth_layers_bmi) {
    if (parameters->soil_storage_model == Layered) {
      stringstream errMsg;
      errMsg << "soil_depth_layers not set in the config file "<< config_file << "\n";
      throw runtime_error(errMsg.str());
    }
  }
  
  if (!is_smcmax_set) {
    stringstream errMsg;
    errMsg << "smcmax not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }
  
  if (!is_b_set) {
    stringstream errMsg;
    errMsg << "b (Clapp-Hornberger's parameter) not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }
  
  if (!is_satpsi_set) {
    stringstream errMsg;
    errMsg << "satpsi not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }

  if (!is_soil_moisture_fraction_depth_set) {
    parameters->soil_moisture_fraction_depth = 0.4; // in meters
  }
  
  if (!is_soil_storage_model_depth_set && parameters->soil_storage_model == Conceptual) {
    stringstream errMsg;
    errMsg << "soil_storage_model_depth not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }

  if (!is_soil_storage_model_set) {
    stringstream errMsg;
    errMsg << "soil_storage_model not set in the config file "<< config_file << ". Options 'conceptual or layered' \n";
    throw runtime_error(errMsg.str());
  }

  
  if (parameters->soil_storage_model == Layered) {
    if (!is_soil_moisture_profile_option_set) {
      stringstream errMsg;
      errMsg << "soil_moisture_profile_option_set key is not set in the config file "<< config_file << ", options = constant or linear \n";
      throw runtime_error(errMsg.str());
    }

    if (!is_water_table_depth_set) {
      parameters->water_table_depth = 6.0;
    }
    
    assert (parameters->num_wetting_fronts > 0);
    assert (parameters->water_table_depth >= 0);
  }

  // check to ensure that options for the topmodel based watertable provided are correct
  if (parameters->soil_storage_model == Topmodel) {
    
    if (is_water_table_based_method_set) {
      if (parameters->water_table_based_method != Flux_based && parameters->water_table_based_method != Deficit_based) {
	stringstream errMsg;
	errMsg << "water_table_based_method key is set in the config file "<< config_file << " with wrong options. Options = flux_based or deficit_based \n";
	throw runtime_error(errMsg.str());
      }
    }
    else {
      stringstream errMsg;
      errMsg << "soil_storage_model is set to Topmodel in the config file "<< config_file << ", but water_table_based_method is not set. Options = iterative or deficit \n";
      throw runtime_error(errMsg.str());
    }
  }
  
  assert (parameters->ncells > 0);

}


void soil_moisture_profile::
SoilMoistureProfileUpdate(struct soil_profile_parameters* parameters)
{
  // update soil moisture fraction
  double thickness = 0.0;
  double soil_moisture_fraction_depth = parameters->soil_moisture_fraction_depth;
  parameters->soil_moisture_fraction = 0.0; // reset to 0 at each timestep
  
  
  if (parameters->soil_storage_model == Conceptual) {
    SoilMoistureProfileFromConceptualReservoir(parameters);
  }
  else if (parameters->soil_storage_model == Layered) {
    SoilMoistureProfileFromLayeredReservoir(parameters);
  }
  else if (parameters->soil_storage_model == Topmodel) {
    SoilMoistureProfileFromWaterTableDepth(parameters);    
  }
  else {
    stringstream errMsg;
    errMsg << "Soil moisture profile OPTION provided in the config file is " << parameters->soil_storage_model
	   << ", valid options are concepttual, layered, and topmodel " <<"\n";
    throw runtime_error(errMsg.str());
  }

  // update the soil storage for the top 2 m to be consistent with NWM; needed only when coupled to Topmodel
  if (parameters->soil_storage_model == Topmodel || parameters->soil_storage_model == Layered) {
    double storage_temp = 0.0;  // temporary storage
    for (int i=0; i<parameters->ncells; i++) {
      if (parameters->soil_z[i] <= parameters->soil_depth_NWM && i == 0) {      
	storage_temp = parameters->soil_moisture_profile[i] * parameters->soil_z[i];
      }
      else if (parameters->soil_z[i] <= parameters->soil_depth_NWM) {
	thickness = parameters->soil_z[i] - parameters->soil_z[i-1];
	storage_temp += parameters->soil_moisture_profile[i] * thickness;
      }
      else {
	break;
      }
      
    }
    parameters->soil_storage = storage_temp;

  }
  
  // compute soil moisture fraction, moisture in the top 40 cm over soil storage in 2 m 
  for (int i=0; i<parameters->ncells; i++) {
    if (parameters->soil_z[i] <= soil_moisture_fraction_depth && i == 0) {      
      parameters->soil_moisture_fraction += parameters->soil_moisture_profile[i] * parameters->soil_z[i];
    }
    else if (parameters->soil_z[i] <= soil_moisture_fraction_depth) {
      thickness = parameters->soil_z[i] - parameters->soil_z[i-1];
      parameters->soil_moisture_fraction += parameters->soil_moisture_profile[i] * thickness;
    }
    else {
      break;
    }
    
  }
  
  parameters->soil_moisture_fraction = fmin(parameters->soil_moisture_fraction, parameters->soil_storage);
  
  if (parameters->soil_storage > 0.0)
    parameters->soil_moisture_fraction /= parameters->soil_storage;
  else
    parameters->soil_moisture_fraction = 0.0;
  
}

/*
  Computes 1D soil moisture profile for conceptual reservoir using Newton-Raphson iterative method
  For detailed decription of the model implemented here, please see README.md on the github repo
  local_variables:
  @param lam  [-] : 1/b (b: pore size distribution)
  @param satpsi_cm [cm] : saturated moisture potential
  @param soil_storage_model_depth [m] : depth of the soil reservoir model (e.g., CFE)
  @param zb  [cm] : bottom of the computational domain
  @param z0  [cm] : bottom of the fictitious domain (to track fictitious water table location)
  @param z1  [cm] : depth of the water table at the previous timestep (known)
  @param z2  [cm] : depth of the water table at the current timestep (unknown, will use the Newton-Raphson method to find it)
  @param soil_storage_max [cm] : maximum soil storage
  @param soil_storage_previous_timestepcm [cm] : soil storage at the previous timestep computed
                                                 locally from the watertable location at the previous timestep
  @param soil_storage_current_timestepcm [cm] : soil storage at the current timestep
  @param tolerance             [cm] : Error tolerance for finding the new root (i.e., water_table_depth)
  @param soil_moisture_profile [-]  : OUTPUT (soil moisture content vertical profile [-])
  ** NOTE: the module needs to be fixed if the CFE and SFT soil depths are different
*/

void soil_moisture_profile::
SoilMoistureProfileFromConceptualReservoir(struct soil_profile_parameters* parameters)
{
  std::string verbosity = parameters->verbosity;
  // converting variables to cm for numerical reasons only
  double satpsi_cm   = parameters->satpsi * 100.;
  double model_depth = parameters->soil_storage_model_depth * 100.;
  double zb          = 0.0;  // bottom of the computational domain
  double z0          = 0.0;  // bottom of the fictitious domain (to track fictitious water table location)
  double zi          = 0.01; // initial guess for the water table location, use Newton-Raphson to find new zi
  double lam         = 1.0/parameters->b;
  double beta        = 1.0 - lam;
  double alpha       = pow(satpsi_cm,lam)/beta; // a constant term obtained in the integration of the soil moisture function
  double tolerance   = 1.0E-6;

  double soil_storage_max = model_depth * parameters->smcmax[0];
  
  double soil_storage_change_per_timestep_cm = fabs(parameters->soil_storage_change_per_timestep * 100.0);
  double soil_storage_current_timestep_cm    = 100.0 * parameters->soil_storage;  // storage at the current timestep
  
  assert(parameters->soil_storage > 0.0); /* to ensure that soil storage is non-zero due to unexpected
					      bugs (either in the models or calibration tools) */
  
  int count = 0;

  /* compute a new profile only if sufficient amount of water is added at this timestep.
     1.0E-4 corresponds to 0.001 mm of water */
  
  if (soil_storage_change_per_timestep_cm > 0.0001 || parameters->init_profile) {
    
    // turn off the flag for times t > 0 
    parameters->init_profile = false;
    
    // check if the storage is greater than the maximum soil storage. if yes, set it to the maximum storage
    if(soil_storage_current_timestep_cm >= soil_storage_max) {
      for(int j=0;j<parameters->ncells;j++)
	parameters->soil_moisture_profile[j] = parameters->smcmax[0];

      parameters->water_table_depth = 0.0;
      return;
    }

    // this won't happen, hopefully, if it does then assign a small non-zero number to the soil moisture profile
    // also to avoid computing water table and over-estimated soil moisture for extremely dry soils
    if(soil_storage_current_timestep_cm <= 1.0E-6) {
      for(int j=0;j<parameters->ncells;j++)
	parameters->soil_moisture_profile[j] = 1.0E-6;

      parameters->water_table_depth = 1000.0; // fictitious water table
      return;
    }
      
    double diff=1000.0; // guess for the initial differnce between the roots
    double f, zi_new, df_dzi;

    // note for stability reasons, the loop is terminated when the water-table depth exceeds 1000m,
    // that means the soil is super dry
    do {
      count++;
      
      // function representing the total amount of soil moisture. 2nd term is the integral of
      // the Clap-Hornberger function (area under the soil moisture curve)

      // fis is integrates the function from z0 to the surface: 1st part: saturated soil between
      // zi and z0; 2nd part: capillary fringe; 3rd: enclosed area between satpis and the surface
      // fib is non-zero for zi <0, zero otherwise. this is the volume of water that needs to be
      // subtracted from "fis" to get the water storage in the computational domain (say top 2 m if soil column has depth 2m)
      double fis = parameters->smcmax[0] * (zi - z0) + parameters->smcmax[0] * satpsi_cm
	           + alpha * parameters->smcmax[0] * ( pow((model_depth-zi),beta) - pow(satpsi_cm,beta) );
      double fib = parameters->smcmax[0] * (zi - z0) + parameters->smcmax[0] * satpsi_cm
	           + alpha * parameters->smcmax[0] * ( pow(abs(zb-zi),beta) - pow(satpsi_cm,beta) );
	

      fib = zi >= 0.0 ? 0.0 : fib;

      f = fis - fib - soil_storage_current_timestep_cm;

      // derivative of f w.r.t zi
      double dfis = parameters->smcmax[0] - alpha * parameters->smcmax[0] * beta * pow((model_depth-zi),beta-1.0);
      double dfib = parameters->smcmax[0] - alpha * parameters->smcmax[0] * beta * pow(fabs(zb-zi),beta-1.0);

      dfib = zi >= 0.0 ? 0.0 : dfib;
      df_dzi = dfis - dfib;
      
      // Newton-Raphson method
      zi_new = zi - f / fmax(df_dzi,1.e-6); // to avoid division by zero
      
      diff=zi_new-zi; // difference betweent the previous and new root
      
      zi=zi_new;     // update the previous root

      z0 = zi >= 0.0 ? zb : zi - satpsi_cm;

      double zi_m = zi/100.; // zi in meters
      // if the water gets below 1000 m, that would mean the soil is super dry and the algorithm
      // may fail to converge in reasonable number of timesteps
      if (zi_m < -1000)
	break;
      
    } while (fabs(diff) > tolerance);

    // water table thickness can be negative and that would be depth of the water table below the
    // depth of the computational domain; probably a better name would be water_table_location
    parameters->water_table_depth = (model_depth - zi)/100.;

    /*******************************************************************/
    // compute the soil moisture profile for the given soil discretization

    for (int i=0; i<parameters->ncells; i++) {
      double z_temp = parameters->water_table_depth - parameters->soil_z[i];
      double theta;
      
      if (parameters->water_table_depth <= parameters->soil_z[i])
	theta = parameters->smcmax[0];
      else
	theta = parameters->smcmax[0] * pow((parameters->satpsi/z_temp),lam);
      
      parameters->soil_moisture_profile[i] = fmin(theta, parameters->smcmax[0]);
      
    }
    
  }

  if (verbosity.compare("high") == 0) {
    std::cout<<"Number of iterations  = "<< count <<"\nWater table depth (m) = "<< parameters->water_table_depth <<"\n";
    PrintSoilMoistureProfile(parameters);
    
    // check compute water in the model domaian
    double total_water = parameters->soil_moisture_profile[0] * parameters->soil_z[0];

    for (int i=1; i<parameters->ncells; i++) {
      total_water += parameters->soil_moisture_profile[i] * (parameters->soil_z[i] - parameters->soil_z[i-1]);
    }
    std::cout<<"Given soil water = "<<parameters->soil_storage<<", Computed soil water = "<<total_water<<"\n";
  }

  for (int i=1; i<parameters->ncells; i++) {
    assert (parameters->soil_moisture_profile[i] <= parameters->smcmax[0]);
    assert (parameters->soil_moisture_profile[i] > 0.0);
  }
}



/*
  - Computes 1D soil moisture profile for layered/calculated reservoir. That is, take layered reservior
    soil moisture and distribute it vertically to the desired soil discretization
  - Two strategies are implemented
  - Constant strategy: A simple technique to map layered values to grids in the soil discretization.
                       That is, all grid cells in the discretization has a constant value within a layer.
		       Note cells at the interface take average value of the layers 
  - Linear strategy: A linear interpolation tehcnique is used to map layered values to grids in the
                     soil discretization. That is, grid cells in the discretization take linearly interpolated
		     value between consecutive layers.
  - Note: the water table location is the thickness of the water table plus saturated capillary head (satpsi)
  - local_variables:
  @param lam                [-]  : 1/b (b: pore size distribution)
  @param soil_depth         [cm] : depth of the soil column
  @param last_layer_depth   [cm] : depth of the last layer from the surface
  @param tolerance          [-]  : tolerance to find location of the water table
*/

void soil_moisture_profile::
SoilMoistureProfileFromLayeredReservoir(struct soil_profile_parameters* parameters)
{
  std::string verbosity = parameters->verbosity;
  int num_wf            = parameters->num_wetting_fronts; //number of wetting fronts
  int num_layers        = parameters->num_layers;

  if (verbosity.compare("high") == 0) {
    std::cerr<<"SoilMoistureProfile: number of wetting fronts = "<<num_wf<<"\n";
    for (int i =0; i <num_wf; i++)
      std::cerr<<"SoilMoistureProfile (input): (depth, water_content) = "<<parameters->soil_depth_wetting_fronts[i]
	       <<", "<<parameters->soil_moisture_wetting_fronts[i]<<"\n";
  }

  parameters->last_layer_depth = parameters->soil_depth_wetting_fronts[num_wf-1];
  
  double lam = 1.0/parameters->b; // pore distribution index
  
  vector<double> z_layers_n(1,0.0);

  
  for (int i=0; i < num_wf; i++)
    z_layers_n.push_back(parameters->soil_depth_wetting_fronts[i]);

  vector<double> smc_column;
  int c = 0;

  // call to find water table function for layered reservoirs
  FindWaterTableLayeredReservoir(parameters);
  
  // piece-wise constant (vertically)
  if (parameters->soil_moisture_profile_option == Constant) {
    
    // loop over all the cells in the discretized column
    for (int i=0; i < parameters->ncells; i++) {

      if (parameters->soil_z[i] <= parameters->last_layer_depth) {

	if (i == 0 && parameters->soil_z[i] > parameters->soil_depth_wetting_fronts[c]) {
	  parameters->soil_moisture_profile[i] = 0.5 * (parameters->soil_moisture_wetting_fronts[c]
							+ parameters->soil_moisture_wetting_fronts[c+1]);
	  c++;
	}
	else if (parameters->soil_z[i] <= parameters->soil_depth_wetting_fronts[c]) {  // cell completely lie within a layer
	  parameters->soil_moisture_profile[i] = parameters->soil_moisture_wetting_fronts[c];
	}
	else { // cell at the interface of layers, so take the mean
	  
	  if (parameters->soil_z[i-1] == parameters->soil_depth_wetting_fronts[c] &&
	      parameters->soil_z[i] <= parameters->soil_depth_wetting_fronts[c+1]) {
	    parameters->soil_moisture_profile[i] = parameters->soil_moisture_wetting_fronts[c+1];
	  }
	  else {
	    parameters->soil_moisture_profile[i] = 0.5 * (parameters->soil_moisture_wetting_fronts[c]
							  + parameters->soil_moisture_wetting_fronts[c+1]);
	  }
	  c++;
	  
	}
      }
      else {
	// extend the profile below the last layer depth, covering some corner cases
	double theta;
	if (parameters->water_table_depth <= parameters->soil_z[i])
	  theta = parameters->smcmax[num_layers-1];
	else {
	  double z_temp = parameters->water_table_depth - parameters->soil_z[i];
	  theta = parameters->smcmax[num_layers-1] * pow((parameters->satpsi/z_temp),lam);
	}
	
	parameters->soil_moisture_profile[i] = theta;
      }
    }

  }
  else if (parameters->soil_moisture_profile_option == Linear ) {
    double value       = 0.0;
    int c = 0; // the parameter c keeps track of wetting fronts
    
    for (int i=0; i < parameters->ncells; i++) {

      if (parameters->soil_z[i] <= parameters->soil_depth_wetting_fronts[0]) {
	parameters->soil_moisture_profile[i] = parameters->soil_moisture_wetting_fronts[0];
	
	if (parameters->soil_z[i+1] > parameters->soil_depth_wetting_fronts[0])
	  c++;
      }
      // linear interpolation between consecutive layers
      else if (parameters->soil_z[i] <= parameters->soil_depth_wetting_fronts[num_wf-1]) {

	if (c == num_wf)
	  break;
	
	value = LinearInterpolation(parameters->soil_z[i], parameters->soil_depth_wetting_fronts[c-1],
				  parameters->soil_depth_wetting_fronts[c], parameters->soil_moisture_wetting_fronts[c-1],
				  parameters->soil_moisture_wetting_fronts[c]);
	
	parameters->soil_moisture_profile[i] = value;

	// if true, we reached the domain depth of the soil moisture profile
	if (i == parameters->ncells-1)
	  break;
	
	if (parameters->soil_z[i+1] > parameters->soil_depth_wetting_fronts[c])
	  c++;
      }
      else {
	// extend the profile below the last layer depth, covering some corner cases
	double theta;
	if (parameters->water_table_depth <= parameters->soil_z[i])
	  theta = parameters->smcmax[num_layers-1];
	else {
	  double z_temp = parameters->water_table_depth - parameters->soil_z[i];
	  theta = parameters->smcmax[num_layers-1] * pow((parameters->satpsi/z_temp),lam);
	}
	
	parameters->soil_moisture_profile[i] = theta;
	
      }
    }
    
  }
  
  if (verbosity.compare("high") == 0) {
    std::cout<<"Water table depth (m) = "<< parameters->water_table_depth <<"\n";
    PrintSoilMoistureProfile(parameters);
  }
  
}

void soil_moisture_profile::
FindWaterTableLayeredReservoir(struct soil_profile_parameters* parameters)
{
  std::string verbosity = parameters->verbosity;
  int num_wf            = parameters->num_wetting_fronts; //number of wetting fronts
  int num_layers        = parameters->num_layers;
  double tolerance      = 1.0e-3;
  double lam = 1.0/parameters->b; // pore distribution index
  bool is_water_table_found = false;
  
  // find and update water table location if it is within the model domain
  if (fabs(parameters->soil_moisture_wetting_fronts[num_wf-1] - parameters->smcmax[num_layers-1]) < tolerance) {
    is_water_table_found = true;
    int j = num_wf-1;
    double z1;
    
    for (int c=num_layers-1; c>=0; c--) {
      z1 = 0.0;
      if ( c != 0)
	z1 = parameters->soil_depth_layers[c-1];
      
      //      std::cerr<<"layer = "<<c<<" "<<parameters->soil_depth_wetting_fronts[j]<<" , "<<parameters->soil_depth_layers[c]<<" "<<z1<<"\n";
      
      while (parameters->soil_depth_wetting_fronts[j] <= parameters->soil_depth_layers[c] && parameters->soil_depth_wetting_fronts[j] > z1) {
	bool is_wf_saturated = fabs(parameters->soil_moisture_wetting_fronts[j] - parameters->smcmax[c]) < tolerance;
	//std::cerr<<"Vx = "<<j<<" "<<is_wf_saturated<<" "<<parameters->soil_depth_wetting_fronts[j]<<" "
	//	 <<parameters->soil_moisture_wetting_fronts[j]<<" "<<parameters->smcmax[c]<<"\n";

	if (is_wf_saturated && j == 0) {
	  parameters->water_table_depth = 0.0;
	  break;
	}
	else if (is_wf_saturated && j > 0) {
	  parameters->water_table_depth = parameters->soil_depth_wetting_fronts[j-1];
	  j--;
	}
	else if (j==0 || !is_wf_saturated)
	  break;
	
      }
    }
  }
  
  // If the watertable is not within the model domain (i.e., the soil is unsaturated in the column),
  // then we find watertable depth as below

  if (!is_water_table_found) {
    //    double target_theta = parameters->soil_moisture_profile[parameters->ncells-1];
    double target_theta = parameters->soil_moisture_wetting_fronts[num_wf-1];
    double dz = 0.00001;
    double theta = 0.0;
    double initial_head = parameters->satpsi; // fully saturated soil
    while (fabs(theta - target_theta) > tolerance) {
      
      // extend the profile below the depth of the last layer
      double z_head = initial_head + dz;
      theta = pow((parameters->satpsi/z_head),lam) * parameters->smcmax[num_layers-1];
      
      assert (theta <= parameters->smcmax[num_layers-1]);
      
      if (theta <= target_theta)
	break;

      if (dz >= 100.0)
	break;

       dz += 0.05;
    }
    
    // watertable depth = domain_depth + depth_to_WT_from_domain_depth + capillary_fringe
    //parameters->water_table_depth = parameters->soil_z[parameters->ncells-1] + dz + parameters->satpsi;
    parameters->water_table_depth = parameters->soil_depth_wetting_fronts[num_wf-1] + dz + parameters->satpsi;
  }

}
/*
void soil_moisture_profile::
ExtendedProfileBelowDomainDepth(struct soil_profile_parameters* parameters) {
  double target_theta = parameters->soil_moisture_profile[parameters->ncells-1];
  double dz           = 0.00001;
  double theta        = 0.0;
  double initial_head = parameters->satpsi; // fully saturated soil
  double tolerance    = 1.0e-3;
  
  while (fabs(theta - target_theta) > tolerance) {
    
    // extend the profile below the depth of the last layer
    double z_head = initial_head + dz;
    theta = pow((parameters->satpsi/z_head),lam) * parameters->smcmax[num_layers-1];
    
    assert (theta <= parameters->smcmax[num_layers-1]);

    if (theta <= target_theta)
      break;

    if (dz >= 100.0)
      break;
    
    dz += 0.05;
  }

}*/

/*
  Computes 1D soil moisture profile for models (e.g. Topmodel) using the water table depth
  For detailed decription of the model implemented here, please see README.md on the github repo
  local_variables:
  @param lam  [-] : 1/b (b: pore size distribution)
  @param satpsi_cm [cm] : saturated moisture potential
  @param soil_moisture_profile [-] : OUTPUT (soil moisture content vertical profile [-])
  @param dt                    [h] : topmodel's timestep
  @param model_depth           [m] : reference depth (datum) for the topmodel

  Note the two methods implemented here are based on
  Method 1 : Eq. (15) in Franchini et al. (1996))
  Method 2 : Eq. (2) in Blazkova et al. (2002)
*/

void soil_moisture_profile::
SoilMoistureProfileFromWaterTableDepth(struct soil_profile_parameters* parameters)
{
  // converting variables to cm for numerical reasons only
  double satpsi_cm   = parameters->satpsi * 100.;
  double lam         = 1.0/parameters->b;
  double model_depth = 600;
  double dt = 1.0;
  double to_cm = 100;
  double theta_fc = parameters->smcmax[0] / 3.0;
  double delta_theta = (parameters->smcmax[0] - theta_fc);
  std::string verbosity = parameters->verbosity;
  
  if (parameters->water_table_based_method == Deficit_based || parameters->init_profile) {
    parameters->water_table_depth = parameters->global_deficit/delta_theta * to_cm + satpsi_cm; // add saturated head to account for capillary fringe
    
    parameters->init_profile = false;  // turn off the flag for times t > 0 
    }
  else if (parameters->water_table_based_method == Flux_based) {
    // Eq. (15) in M. Franchini et al. Journal of Hydrology 175 (1996) 293-338
    parameters->water_table_depth -= (parameters->Qv_topmodel - parameters->Qb_topmodel)/parameters->cat_area * dt * to_cm;
  }
  
  // check if the storage is greater than the maximum soil storage. if yes, set it to the maximum storage
  if(parameters->water_table_depth == 0.0) {
    for(int j=0;j<parameters->ncells;j++)
      parameters->soil_moisture_profile[j] = parameters->smcmax[0];
    
    return;
  }
  
  /*******************************************************************/
  // get a high resolution moisture profile that will be mapped on the desired soil discretization
    
  int z_hres = 1000;
  double *smct_temp = new double[z_hres];
  double *z_temp = new double[z_hres];

  double zi = model_depth - parameters->water_table_depth; // thickness of the water table depth (bottom to top)
  
  // we have the new water table location now, so let's compute the soil moisture curve now
  double z = zi + satpsi_cm;
  double dz_v = (model_depth - zi - satpsi_cm)/z_hres; // vertical spacing over the depth (from zi+satpsi to the surface)
  double z_head = satpsi_cm;  // capillary head (psi) in the soil moisture function

  for (int i=0;i<z_hres;i++) {
    z_head += dz_v;
    smct_temp[i] = parameters->smcmax[0] * pow((satpsi_cm/z_head),lam) ;

    z+=dz_v;
    z_temp[i] = z;
  }
    
  // map the high resolution soil moisture curve to the soil discretization depth that is provided in the config file
  for (int i=0; i<parameters->ncells; i++) {
    for (int j=0; j<z_hres; j++) {
      if ( (model_depth - parameters->soil_z[i]*100) <= (zi + satpsi_cm) ) {
	parameters->soil_moisture_profile[i] = parameters->smcmax[0];
      }
      else if (z_temp[j]  >= (model_depth - parameters->soil_z[i]*100) ) {
	parameters->soil_moisture_profile[i] = smct_temp[j];
	break;
      }
    }
  }
  
  
  if (verbosity.compare("high") == 0) {
    std::cout<<"Water table depth (m) = "<< parameters->water_table_depth <<"\n";
    PrintSoilMoistureProfile(parameters);
  }

}


double soil_moisture_profile::
LinearInterpolation(double z, double z1, double z2, double t1, double t2)
{
  double m = (t2 - t1) / (z2 - z1);
  return t1 + m * (z - z1);
}


/*
Reads 1D data from the config file
- used for reading soil discretization (1D)
- used for reading layers depth from the surface if model `layered` is chosen
*/

vector<double> soil_moisture_profile::
ReadVectorData(string param_name, string param_value)
{
  int pos =0;
  string delimiter = ",";
  vector<double> values(0.0);
  string z1 = param_value;
  
  if (z1.find(delimiter) == string::npos) {
    double v = stod(z1);
    if (v <= 0.0) {
      stringstream errMsg;
      errMsg << "Input provided in the config file for parameter "<< param_name << " is " << v << ". It should be positive."<< "\n";
      throw runtime_error(errMsg.str());
    }
    
    values.push_back(v);
    
  }
  else {
    while (z1.find(delimiter) != string::npos) {
      pos = z1.find(delimiter);
      string z_v = z1.substr(0, pos);

      values.push_back(stod(z_v.c_str()));
      
      z1.erase(0, pos + delimiter.length());
      if (z1.find(delimiter) == string::npos)
	values.push_back(stod(z1));
    }
  }
  
  return values;
}

void soil_moisture_profile::
PrintSoilMoistureProfile(struct soil_profile_parameters* parameters)
{
  for (int i=0; i<parameters->ncells; i++)
    std::cout<<"soil_moisture (z, value) = "<< parameters->soil_z[i]<<", "<<parameters->soil_moisture_profile[i]<<"\n";
}
  

#endif
