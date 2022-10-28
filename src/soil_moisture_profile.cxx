#ifndef SMP_C_INCLUDED
#define SMP_C_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cmath> // fmin and fmax (C)
//#include <algorithm>  // for std::min and std::max (C++)

#include "../include/soil_moisture_profile.hxx"

enum {Conceptual=1, Layered=2};
enum {Constant=1, Linear=2};

// soil_moisture_profile is the namespacing


void soil_moisture_profile::
SoilMoistureProfile(string config_file, struct soil_profile_parameters* parameters)
{
  
  InitFromConfigFile(config_file, parameters);

  if (parameters->soil_storage_model == Conceptual) {
    parameters->shape[0] = parameters->ncells;
    parameters->shape[1] = 1;
  }
  else if (parameters->soil_storage_model == Layered) {
    parameters->shape[0] = 1;
    parameters->shape[1] = parameters->max_ncells_layered; // note this will be set dynamically at each timestep
    
  }

  parameters->soil_moisture_profile = new double[parameters->ncells];
  
  // the following two will be reallocated at each time step, if coupled
  if (parameters->soil_depths_layered_bmi) {
    parameters->soil_moisture_layered = new double[parameters->shape[1]];
    parameters->soil_depths_layered = new double[parameters->shape[1]];
  }
  
  parameters->shape[2] = 1;
  parameters->spacing[0] = 1.;
  parameters->spacing[1] = 1.;
  parameters->origin[0] = 0.;
  parameters->origin[1] = 0.;
  parameters->soil_storage = 0.0;
  parameters->soil_storage_change_per_timestep = 0.0;
  parameters->init_profile = true;
  parameters->ncells_layered = 1;
}

/*
Read and initialize values from configuration file
@input - soil_z   (1D)  :   soil discretization; array of depths from the surface [m]
@input - layers_z  (1D) : depth of each layer from the surface [m]
@input - bb  (double)  : pore size distribution [-], beta exponent on Clapp-Hornberger (1978)
@input - satpsi  (double) : saturated capillary head (saturated moisture potential) [m]
@input - ncells  (int) : number of cells of the discretized soil column
@input - ncells_layered : number of layers (wetting fronts) for the layered model
@input - max_ncells_layered : maximum number of layers (wetting fronts)
@input - soil_storage_model (string) : Conceptual or Layered soil reservoir models
@input - soil_moisture_profile_option (string) : valid only when layered soil reservoir model is chosen; option include `constant` or `linear`. The option `constant` assigns a constant value to discretized cells within each layer, `linear` option linearly interpolate values between layers and interpolated values are assigned to the soil discretization
@params - input_var_names_model (1D) : dynamically sets model inputs to be used in the bmi input_var_names
@input water_table_depth (double)    : water table depth, default is 6 m for layered model, conceptual model computes its own
*/


void soil_moisture_profile::
InitFromConfigFile(string config_file, struct soil_profile_parameters* parameters)
{ 
  ifstream fp;
  fp.open(config_file);
  
  bool is_soil_z_set = false;
  bool is_soil_depths_layered_set = false;
  bool is_smcmax_set = false;
  bool is_bb_set = false;
  bool is_satpsi_set = false;
  bool is_soil_storage_model_set = false;
  bool is_soil_moisture_layered_option_set = false; // option for linear or piece-wise constant layered profile
  bool is_soil_storage_model_depth_set = false;
  bool is_max_ncells_layered_set = false;
  bool is_water_table_depth_set = false;
  
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
      vector<double> vec = ReadVectorData(param_value);
      
      parameters->soil_z = new double[vec.size()];
      
      for (unsigned int i=0; i < vec.size(); i++)
      	parameters->soil_z[i] = vec[i];
      
      parameters->ncells = vec.size();
      parameters->soil_depth = parameters->soil_z[parameters->ncells-1];
      is_soil_z_set = true;
      continue;
    }
    else if (param_key == "soil_depths_layered") {
      if (param_value == "bmi" || param_value == "BMI") {
	parameters->soil_depths_layered_bmi = true;
      }
      else {
	vector<double> vec = ReadVectorData(param_value);
	parameters->soil_depths_layered = new double[vec.size()];
	
	for (unsigned int i=0; i < vec.size(); i++)
	  parameters->soil_depths_layered[i] = vec[i];
	
	parameters->ncells_layered = vec.size();
	parameters->last_layer_depth = parameters->soil_depths_layered[parameters->ncells_layered-1];
	is_soil_depths_layered_set = true;
	parameters->soil_depths_layered_bmi = false;
      }
      continue;
    }
    else if (param_key == "soil_params.smcmax") {
      parameters->smcmax = stod(param_value);
      is_smcmax_set = true;
      continue;
    }
    else if (param_key == "soil_params.b") {
      parameters->bb = stod(param_value);
      assert (parameters->bb > 0);
      is_bb_set = true;
      continue;
    }
    else if (param_key == "soil_params.satpsi") {
      parameters->satpsi = stod(param_value);
      is_satpsi_set = true;
      continue;
    }
    else if (param_key == "soil_storage_model") {
      if ( param_value == "Conceptual" || param_value == "conceptual")
	parameters->soil_storage_model = Conceptual;
      else if (param_value == "layered" || param_value == "Layered") 
	parameters->soil_storage_model = Layered;

      is_soil_storage_model_set = true;
      continue;
    }
    else if (param_key == "soil_moisture_layered_option") {  //Soil moisture profile option
      if (param_value == "Constant" || param_value == "constant")
	parameters->soil_moisture_layered_option = Constant;
      else if (param_value == "Linear" || param_value == "linear")
	parameters->soil_moisture_layered_option = Linear;
      is_soil_moisture_layered_option_set = true;
      continue;
    }
    else if (param_key == "soil_storage_depth") {
      parameters->soil_storage_model_depth = stod(param_value);
      assert (parameters->soil_storage_model_depth > 0);
      is_soil_storage_model_depth_set = true;
      continue;
    }
    else if (param_key == "max_num_cells_layered") {
      parameters->max_ncells_layered =  stod(param_value);
      assert (parameters->max_ncells_layered > 0);
      is_max_ncells_layered_set = true;
      continue;
    }
    else if (param_key == "water_table_depth") {
      parameters->water_table_depth = stod(param_value);
      is_water_table_depth_set = true;
      continue;
    }
    else if (param_key == "verbosity") {
      verbosity = param_value;
      continue;
    }
  }
  
  fp.close();
  
  if (!is_soil_z_set) {
    stringstream errMsg;
    errMsg << "soil_z not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }

  
  if (!is_soil_depths_layered_set && !parameters->soil_depths_layered_bmi) {
    if (parameters->soil_storage_model == Layered) {
      stringstream errMsg;
      errMsg << "soil_depths_layered not set in the config file "<< config_file << "\n";
      throw runtime_error(errMsg.str());
    }
  }

  if (!is_smcmax_set) {
    stringstream errMsg;
    errMsg << "smcmax not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }
  
  if (!is_bb_set) {
    stringstream errMsg;
    errMsg << "bb (Clapp-Hornberger's parameter) not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }
  
  if (!is_satpsi_set) {
    stringstream errMsg;
    errMsg << "satpsi not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
    }
  
  if (!is_soil_storage_model_depth_set && parameters->soil_storage_model == Conceptual) {
    stringstream errMsg;
    errMsg << "soil_storage_model_depth not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }

  
  if (parameters->soil_storage_model == Layered) {
    if (!is_soil_moisture_layered_option_set) {
      stringstream errMsg;
      errMsg << "soil_moisture_layered_option_set key is not set in the config file "<< config_file << ", options = constant or linear \n";
      throw runtime_error(errMsg.str());
    }
    
    if (!is_max_ncells_layered_set) {
      parameters->max_ncells_layered = 30;
    }

    if (!is_water_table_depth_set) {
      parameters->water_table_depth = 6.0;
    }
    
    assert (parameters->max_ncells_layered > 0);
    assert (parameters->water_table_depth >= 0);
  }

  assert (parameters->ncells > 0);

}


void soil_moisture_profile::
SoilMoistureProfileUpdate(struct soil_profile_parameters* parameters)
{
  if (parameters->soil_storage_model == Conceptual) {
    SoilMoistureProfileFromConceptualReservoir(parameters);
  }
  else if (parameters->soil_storage_model == Layered) {
    SoilMoistureProfileFromLayeredReservoir(parameters);
  }
  else {
    stringstream errMsg;
    errMsg << "Soil moisture profile OPTION provided in the config file is " << parameters->soil_storage_model<< ", which should be either \'concepttual\' or \'layered\' " <<"\n";
    throw runtime_error(errMsg.str());
  }

}

/*
  Computes 1D soil moisture profile for conceptual reservoir using Newton-Raphson iterative method
  For detailed decription of the model implemented here, please see README.md on the github repo
  local_variables:
  @param lam  [-] : 1/bb (bb: pore size distribution)
  @param satpsi_cm [cm] : saturated moisture potential
  @param soil_storage_model_depth [m] : depth of the soil reservoir model (e.g., CFE)
  @param zb  [cm] : bottom of the computational domain
  @param z0  [cm] : bottom of the fictitious domain (to track fictitious water table location)
  @param z1  [cm] : depth of the water table at the previous timestep (known)
  @param z2  [cm] : depth of the water table at the current timestep (unknown, will use the Newton-Raphson method to find it)
  @param soil_storage_max [cm] : maximum soil storage
  @param soil_storage_previous_timestepcm [cm] : soil storage at the previous timestep computed locally from the watertable location at the previous timestep
  @param soil_storage_current_timestepcm [cm] : soil storage at the current timestep
  @param tolerance             [cm] : Error tolerance for finding the new root (i.e., water_table_depth)
  @param soil_moisture_profile [-]  : OUTPUT (soil moisture content vertical profile [-])
  ** NOTE: the module needs to be fixed if the CFE and SFT soil depths are different
*/

void soil_moisture_profile::
SoilMoistureProfileFromConceptualReservoir(struct soil_profile_parameters* parameters)
{
  // converting variables to cm for numerical reasons only
  double satpsi_cm = parameters->satpsi * 100.;
  double model_depth = parameters->soil_storage_model_depth * 100.;
  double zb = 0.0; // bottom of the computational domain
  double z0 = 0.0; // bottom of the fictitious domain (to track fictitious water table location)
  
  double zi = 0.01; // initial guess for the water table location, use Newton-Raphson to find new zi
  double soil_storage_max = model_depth * parameters->smcmax;
  
  double soil_storage_change_per_timestep_cm = fabs(parameters->soil_storage_change_per_timestep * 100.0);
  double soil_storage_current_timestep_cm = 100.0 * parameters->soil_storage;  // storage at the current timestep
  
  assert(parameters->soil_storage >= 0.0); // to ensure that soil storage is non-negative due to unexpected bugs in cfe (or any other conceptual models)
  
  double lam = 1.0/parameters->bb;
  double beta = 1.0 - lam;
  double alpha = pow(satpsi_cm,lam)/beta; // a constant term obtained in the integration of the soil moisture function
  double tolerance = 0.000001;

  int count = 0;

  // compute a new profile only if sufficient amount of water is added at this timestep. 1.0E-4 corresponds to 0.001 mm of water
  
  if (soil_storage_change_per_timestep_cm > 0.0001 || parameters->init_profile) {
    
    // turn of the flag for times t > 0 
    parameters->init_profile = false;
    
    // check if the storage is greater than the maximum soil storage. if yes, set it to the maximum storage
    if(soil_storage_current_timestep_cm >= soil_storage_max) {
      for(int j=0;j<parameters->ncells;j++)
	parameters->soil_moisture_profile[j] = parameters->smcmax;
      return;
    }
    
    double diff=1000.0; // guess for the initial differnce between the roots
  
    double f, zi_new, df_dzi;
    
    do {
      count++;
      
      if(count>10000) {
	throw runtime_error("No convergence loop count: after 10000 iterations!");
      }
      
      // function representing the total amount of soil moisture. 2nd term is the integral of the Clap-Hornberger function (area under the soil moisture curve)

      // fis is integrates the function from z0 to the surface: 1st part: saturated soil between zi and z0; 2nd part: capillary fringe; 3rd: enclosed area between satpis and the surface
      // fib is non-zero for zi <0, zero otherwise. this is the volume of water that needs to be subtracted from "fis" to get the water storage in the computational domain (say top 2 m if soil column has depth 2m)
      double fis = parameters->smcmax * (zi - z0) + parameters->smcmax * satpsi_cm + alpha * parameters->smcmax * ( pow((model_depth-zi),beta) - pow(satpsi_cm,beta) );
      double fib = parameters->smcmax * (zi - z0) + parameters->smcmax * satpsi_cm + alpha * parameters->smcmax * ( pow(abs(zb-zi),beta) - pow(satpsi_cm,beta) );
	

      fib = zi >= 0.0 ? 0.0 : fib;

      f = fis - fib - soil_storage_current_timestep_cm;

      // derivative of f w.r.t zi
      double dfis = parameters->smcmax - alpha * parameters->smcmax * beta * pow((model_depth-zi),beta-1.0);
      double dfib = parameters->smcmax - alpha * parameters->smcmax * beta * pow(fabs(zb-zi),beta-1.0);

      dfib = zi >= 0.0 ? 0.0 : dfib;
      df_dzi = dfis - dfib;
      
      // Newton-Raphson method
      zi_new = zi - f / fmax(df_dzi,1.e-6); // to avoid division by zero
      
      diff=zi_new-zi; // difference betweent the previous and new root
      
      zi=zi_new;     // update the previous root

      z0 = zi >= 0.0 ? zb : zi - satpsi_cm;

      double zi_m = zi/100.; // zi in meters
      // if the water gets below 1000 m, that would mean the soil is super dry and the algorithm may fail to converge in reasonable number of timesteps
      if (zi_m < -1000)
	break;
      
    } while (fabs(diff) > tolerance);

    // water table thickness can be negative and that would be depth of the water table below the depth of the computational domain; probably a better name would be water_table_location
    parameters->water_table_depth = (model_depth - zi)/100.;
    
    /*******************************************************************/
    // get a high resolution moisture profile that will be mapped on the desired soil discretization
    
    int z_hres = 1000;
    double *smct_temp = new double[z_hres];
    double *z_temp = new double[z_hres];
    
    // we have the new water table location now, so let's compute the soil moisture curve now
    double z = zi + satpsi_cm;
    double dz_v = (model_depth - zi - satpsi_cm)/z_hres; // vertical spacing
    double z_head = satpsi_cm; // input variable to the soil moisture function
    
    for (int i=0;i<z_hres;i++) {
      z_head += dz_v;
      smct_temp[i] = parameters->smcmax * pow((satpsi_cm/z_head),lam) ;
      
      z+=dz_v;
      z_temp[i] = z;
    }
    
    // map the high resolution soil moisture curve to the soil discretization depth that is provided in the config file
    for (int i=0; i<parameters->ncells; i++) {
      for (int j=0; j<z_hres; j++) {
	if ( (model_depth - parameters->soil_z[i]*100) <= (zi + satpsi_cm) )
	  parameters->soil_moisture_profile[i] = parameters->smcmax;
	else if (z_temp[j]  >= (model_depth - parameters->soil_z[i]*100) ) {
	  parameters->soil_moisture_profile[i] = smct_temp[j];
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
  @param tolerance         [-]  : tolerance to find location of the water table
*/

void soil_moisture_profile::
SoilMoistureProfileFromLayeredReservoir(struct soil_profile_parameters* parameters)
{
  int ncells_layered = parameters->ncells_layered; //number of wetting fronts
  double tolerance = 1.0e-4;

  if (verbosity.compare("high") == 0) {
    std::cerr<<"SoilMoistureProfile: number of wetting fronts = "<<ncells_layered<<"\n";
    for (int i =0; i <ncells_layered; i++)
      std::cerr<<"SoilMoistureProfile (input): (depth, water_content) = "<<parameters->soil_depths_layered[i]<<", "<<parameters->soil_moisture_layered[i]<<"\n";
  }

  parameters->last_layer_depth = parameters->soil_depths_layered[ncells_layered-1];
  
  double lam = 1.0/parameters->bb; // pore distribution index
  
  vector<double> z_layers_n(1,0.0);

  
  for (int i=0; i < ncells_layered; i++)
    z_layers_n.push_back(parameters->soil_depths_layered[i]);

  vector<double> smc_column;
  int c = 0;
  double delta = 0.0;

  // piece-wise constant (vertically)
  if (parameters->soil_moisture_layered_option == Constant) {
    bool layers_flag=true;
    
    // loop over all the cells in the discretized column
    for (int i=0; i < parameters->ncells; i++) {

      if (parameters->soil_z[i] <= parameters->last_layer_depth) {

	if (i == 0 && parameters->soil_z[i] > parameters->soil_depths_layered[c]) {
	  parameters->soil_moisture_profile[i] = 0.5*(parameters->soil_moisture_layered[c] + parameters->soil_moisture_layered[c+1]);
	  c++;
	}
	else if (parameters->soil_z[i] <= parameters->soil_depths_layered[c]) {  // cell completely lie within a layer
	  parameters->soil_moisture_profile[i] = parameters->soil_moisture_layered[c];
	}
	else { // cell at the interface of layers, so take the mean
	  
	  if (parameters->soil_z[i-1] == parameters->soil_depths_layered[c] && parameters->soil_z[i] <= parameters->soil_depths_layered[c+1]) {
	    
	  parameters->soil_moisture_profile[i] = parameters->soil_moisture_layered[c+1];
	  }
	  else {
	    parameters->soil_moisture_profile[i] = 0.5*(parameters->soil_moisture_layered[c] + parameters->soil_moisture_layered[c+1]);
	  }
	  c++;
	  
	}
      }
      else { // extend the profile below the last layer depth
	
	//double zz = fmax(parameters->soil_depth - parameters->soil_z[i], 1.0e-4);
	double zz = fmax(parameters->water_table_depth - parameters->soil_z[i], 1.0e-4);
	double theta = delta + pow((parameters->satpsi/zz),lam) * parameters->smcmax;
	
	if (layers_flag) { // Ensure continuity of soil moisture at the interface of the depth of the last layer and profile below that depth of the last layer
	  delta = 0.0 - theta;
	  theta = theta + delta;
	  layers_flag=false;
	}

	double theta1 = parameters->soil_moisture_layered[ncells_layered-1] +  theta;
	
	parameters->soil_moisture_profile[i] = fmin(parameters->smcmax, theta1);

      }
      
    }

    if (verbosity.compare("high") == 0) {
      for (int j=0; j< parameters->ncells; j++)
	cerr<<"SoilMoistureProfile (output): (depth, water_content) = "<<parameters->soil_z[j] <<", "<<parameters->soil_moisture_profile[j]<<"\n";
    }
  }
  else if (parameters->soil_moisture_layered_option == Linear ) {
    
    bool layers_flag=true;
    double t_v=0.0;
   
    parameters->soil_moisture_profile[0] = parameters->soil_moisture_layered[0]; // first cell gets the top-layer soil moisture

    for (int i=1; i < parameters->ncells; i++) {

      // linear interpolation between consecutive layers
      if (parameters->soil_z[i] <= parameters->soil_depths_layered[c] && c < ncells_layered-1) {
	t_v = LinearInterpolation(z_layers_n[c], z_layers_n[c+1], parameters->soil_moisture_layered[c], parameters->soil_moisture_layered[c+1], parameters->soil_z[i]);
	parameters->soil_moisture_profile[i] = t_v;
	
	if (parameters->soil_z[i+1] > parameters->soil_depths_layered[c]) // the parameter c keeps track of the layer
	  c++;
      }
      else {
	// extend the profile below the depth of the last layer
	double zz = fmax(parameters->water_table_depth - parameters->soil_z[i], 1.0e-4);
	double theta = delta + pow((parameters->satpsi/zz),lam)*parameters->smcmax;
	
	if (layers_flag) { // Ensure continuity of soil moisture at the interface of the depth of the last layer and profile below that depth of the last layer
	  delta = 0.0 - theta;
	  theta = theta + delta;
	  layers_flag = false;
	}
	
	double theta1 = parameters->soil_moisture_layered[ncells_layered-1] + theta;

	parameters->soil_moisture_profile[i] = fmin(parameters->smcmax, theta1);
	
      }
    }

    if (verbosity.compare("high") == 0) {
      for (int j=0; j< parameters->ncells; j++)
	cerr<<"SoilMoistureProfile (output): (depth, water_content) = "<<parameters->soil_z[j] <<", "<<parameters->soil_moisture_profile[j]<<"\n";
    }
    
  }
  else {
    stringstream errMsg;
    errMsg << "Soil moisture profile "<< parameters->soil_storage_model << " works with options \'constant\' and \'linear\'. Provide at least one option in the config file "<<"\n";
    throw runtime_error(errMsg.str());
  }

  // for vis comparison with the python version
  //  for (int j=0; j< parameters->ncells; j++)
  //  cout<<parameters->soil_z[j] <<" "<<parameters->soil_moisture_profile[j]<<"\n";
  
  // find and update water table location
  for (int j=0; j< parameters->ncells; j++) {
    if ( fabs(parameters->soil_moisture_profile[j] - parameters->smcmax) < tolerance) {
      parameters->water_table_depth = parameters->soil_z[j];
      break;
    }
  }
  
}


double soil_moisture_profile::
LinearInterpolation(double z1, double z2, double t1, double t2, double z)
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
ReadVectorData(string key)
{
  int pos =0;
  string delimiter = ",";
  vector<double> value(0.0);
  string z1 = key;
  
  if (z1.find(delimiter) == string::npos) {
    double v = stod(z1);
    if (v == 0.0) {
      stringstream errMsg;
      errMsg << "soil_z (depth of soil reservior) should be greater than zero. It it set to "<< v << " in the config file "<< "\n";
      throw runtime_error(errMsg.str());
    }
    
    value.push_back(v);
    
  }
  else {
    while (z1.find(delimiter) != string::npos) {
      pos = z1.find(delimiter);
      string z_v = z1.substr(0, pos);

      value.push_back(stod(z_v.c_str()));
      
      z1.erase(0, pos + delimiter.length());
      if (z1.find(delimiter) == string::npos)
	value.push_back(stod(z1));
    }
  }
  
  return value;
}


/*
soil_moisture_profile::SoilMoistureProfile::
~SoilMoistureProfile()
{}
*/
#endif
