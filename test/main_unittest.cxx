/*
  @author Ahmad Jan (ahmad.jan@noaa.gov)
  - Includes unit test for bmi components and run model for a timestep to compute soil moisture profile and water table lcoation
  - benchmark water table location is compared against computed water table
  - test passed criteria : pass bmi functions, i.e., bmi functions returns/sets expected values/behavior and computed water table matches unittest water table depth
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "../bmi/bmi.hxx"
#include "../include/bmi_soil_moisture_profile.hxx"
#include "../include/soil_moisture_profile.hxx"

#define FAILURE 0
#define VERBOSITY 1

#define GREEN "\033[32m"
#define RED   "\033[31m"
#define BLUE  "\033[34m"
#define RESET "\033[0m"


int main(int argc, char *argv[])
{
  BmiSoilMoistureProfile model, model_layered;

  if (argc != 4) {
    printf("Usage: ./run_unittest.sh CONFIGURATION_FILE1 CONFIGURATION_FILE2 CONFIGURATION_FILE3\n");
    printf("Number of arguments provided = %d \n", argc);
    printf("Runs soil moisture profile model through its BMI with THREE configuration files.\n");
    return FAILURE;
  }

  std::cout<<"\n**************** BEGIN SoilFreezeThaw BMI UNIT TEST *******************\n";
  
  model.Initialize(argv[1]);
  model_layered.Initialize(argv[2]);

  std::cout<<"\n**************** TEST VALUES ************************************\n";
  int nz = 4;
  bool test_status = true;
  int num_input_vars = 8;
  int num_output_vars = 3;
  
  std::vector<string> bmi_input_vars = {"soil_storage", "soil_storage_change", "soil_moisture_wetting_fronts",
					"soil_depth_wetting_fronts", "num_wetting_fronts", "Qb_topmodel", "Qv_topmodel",
					"global_deficit"};
  std::vector<string> bmi_output_vars = {"soil_moisture_profile", "soil_water_table","soil_moisture_fraction"};
  
  int nbytes_input[] = {sizeof(double), sizeof(double), sizeof(double), sizeof(double), sizeof(int),
			sizeof(double), sizeof(double), sizeof(double)};
  int nbytes_output[] = {int(nz * sizeof(double)), sizeof(double), sizeof(double)};
  //double soil_moisture_profile[] = {0.389,0.396,0.397,0.397}; // total_moisture_content
  
  std::cout<<"Num cells:           "<<nz<<"\n";
  std::cout<<"Num input vars:      "<<num_input_vars<<"\n";
  std::cout<<"Num output vars:     "<<num_output_vars<<"\n";

  std::cout<<"\nPulling information from BMI\n************************************\n";
  std::string model_name;
  int count_in = 0;
  int count_out = 0;
  std::vector<std::string> names_in;
  std::vector<std::string> names_out;
  

  // Test get_component_name()
  model_name = model.GetComponentName();
  if (VERBOSITY)
    std::cout<<"Model name: "<< model_name <<"\n";

  // Test GetInputItemCount
  count_in = model.GetInputItemCount();
  if (VERBOSITY)
    std::cout<<"Input item count: "<< count_in<<"\n";
  if (count_in == num_input_vars)
    test_status &= true;
  else {
    test_status &= false;
    std::string passed = test_status == true ? "Yes" : "No";
    std::cout<<"Test passed: "<<passed<<"\n";
    std::stringstream errMsg;
    errMsg << "Number of input variables are different. "<< count_in << " != "<< num_input_vars << "\n";
    throw std::runtime_error(errMsg.str());
  }

  // Test GetInputVarNames 
  names_in = model.GetInputVarNames();
  if (VERBOSITY) {
    std::cout<<"Input variable names \n";
    for (int i=0; i<count_in; i++)
      std::cout<<i<<" "<<names_in[i]<<"\n";
  }

  std::cout<<"**************************************** \n";
  // Test GetOutputItemCount
  count_out = model.GetOutputItemCount();
  if (VERBOSITY)
    std::cout<<"Output item count: "<< count_out<<"\n";

  // Test GetOutputVarNames
  names_out = model.GetOutputVarNames();
  if (VERBOSITY) {
    std::cout<<"Output variable names "<<names_out.size()<<"\n";
    for (int i=0; i<count_out; i++)
      std::cout<<i<<" "<<names_out[i]<<"\n";
  }
  if (count_out == num_output_vars)
    test_status &= true;
  else {
    test_status &= false;
    std::string passed = test_status == true ? "Yes" : "No";
    std::cout<<"Test passed: "<<passed<<"\n";
    std::stringstream errMsg;
    errMsg << "Number of output variables are different. "<< count_out <<" != "<< num_output_vars <<"\n";
    throw std::runtime_error(errMsg.str());
  }
  
  // Test BMI: VARIABLE INFORMATION FUNCTIONS
  std::cout<<"\n**************** TEST BMI VARIABLE INFORMATION FUNCTIONS\n***************************\n";

  int grid, itemsize, nbytes;
  std::string location;
  std::string units;
  std::string vartype;
  
  // Loop through both input and output variables and call GetVar* functions
  for (int i=0; i<count_in; i++) {
    std::string var_name = names_in[i];
    if (VERBOSITY)
      std::cout<<"Input var_name: "<< var_name <<"\n";

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_grid()
    grid = model.GetVarGrid(var_name);
    if (VERBOSITY)
      std::cout<<"Grid: "<< grid <<"\n";

    if (grid >=0)
      test_status &= true;
    else {
      test_status &= false;
      std::string passed = test_status == true ? "Yes" : "No";
      std::cout<<"Test passed: "<<passed<<"\n";
      std::stringstream errMsg;
      errMsg << "grid < 0 \n";
      throw std::runtime_error(errMsg.str());
    }

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_itemsize()
    itemsize = model.GetVarItemsize(var_name);
    if (VERBOSITY)
      std::cout<<"Itemsize: "<< itemsize <<"\n";
    if (itemsize >0)
      test_status &= true;
    else {
      test_status &= false;
      std::string passed = test_status == true ? "Yes" : "No";
      std::cout<<"Test passed: "<<passed<<"\n";
      std::stringstream errMsg;
      errMsg << "itemsize < 0 \n";
      throw std::runtime_error(errMsg.str());
    }

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_location()
    location = model.GetVarLocation(var_name);
    if ( location == "") return FAILURE;
    if (VERBOSITY)
      std::cout<<" location: "<< location<<"\n";
    if (location == "")
      test_status &= false;
    else
      test_status &= true;
	
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_units()
    units = model.GetVarUnits(var_name);
    if (VERBOSITY)
      std::cout<<" units: ["<< units <<"]\n";

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_type()
    vartype = model.GetVarType(var_name);
    if (VERBOSITY)
      std::cout<<" type: "<< vartype <<"\n";
    if (location == "")
      test_status &= false;
    else
      test_status &= true;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // get_var_nbytes()
    nbytes = model.GetVarNbytes(var_name);
    if (nbytes == 0) return FAILURE;
    if (VERBOSITY)
      std::cout<<" nbytes: "<< nbytes <<"\n";

    if (var_name == bmi_input_vars[i]) {
      if (nbytes == nbytes_input[i])
	test_status &= true;
      else {
	test_status &= false;
	std::string passed = test_status == true ? "Yes" : "No";
	std::cout<<"Test passed: "<<passed<<"\n";
	std::stringstream errMsg;
	errMsg << "Number of bytes for input var "<<var_name<< " should be "<<nbytes_input[i]<<" but bmi returned "<<nbytes<<"\n";
	throw std::runtime_error(errMsg.str());
      }
    }
    else {
      std::stringstream errMsg;
      errMsg << "Input variable name "<< var_name <<" not in the bmi_input_vars list. \n";
      throw std::runtime_error(errMsg.str());

    }
  }

  if (VERBOSITY)
    std::cout<<"\n*****************************************\n";
  
  for (int i=0; i<count_out; i++) {
    std::string var_name = names_out[i];
    if (VERBOSITY)
      std::cout<<"Output var_name: "<< var_name <<"\n";
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_grid() 
    grid = model.GetVarGrid(var_name);
    std::cout<<grid<<"\n";
    if (grid == -1) return -1;
    if (VERBOSITY)
      std::cout<<"Grid: "<< grid <<"\n";

    if (grid >=0)
      test_status &= true;
    else
      test_status &= false;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_itemsize()
    itemsize = model.GetVarItemsize(var_name);
    if (itemsize == 0) return FAILURE;
    if (VERBOSITY)
      std::cout<<"Itemsize: "<< itemsize <<"\n";

    if (itemsize >0)
      test_status &= true;
    else
      test_status &= false;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_location()
    location = model.GetVarLocation(var_name);
    if ( location == "") return FAILURE;
    if (VERBOSITY)
      std::cout<<" location:"<< location<<"\n";

    if (location == "")
      test_status &= false;
    else
      test_status &= true;
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_units()
    units = model.GetVarUnits(var_name);
    if (VERBOSITY)
      std::cout<<" units: ["<< units <<"]\n";
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_var_type()
    vartype = model.GetVarType(var_name);
    if (vartype == "") return FAILURE;
    if (VERBOSITY)
      std::cout<<" type: "<< vartype <<"\n";
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // get_var_nbytes()
    nbytes = model.GetVarNbytes(var_name);
    if (nbytes == 0) return FAILURE;
    if (VERBOSITY)
      std::cout<<" nbytes: "<< nbytes<<"\n";

    if (var_name == bmi_output_vars[i]) {
      if (nbytes == nbytes_output[i])
	test_status &= true;
      else {
	test_status &= false;
	std::string passed = test_status == true ? "Yes" : "No";
	std::cout<<"Test passed: "<<passed<<"\n";
	std::stringstream errMsg;
	errMsg << "Number of bytes for output var"<<var_name<< " should be "<<nbytes_output[i]<<"\n";
	throw std::runtime_error(errMsg.str());
      }
    }
    else {
      std::stringstream errMsg;
      errMsg << "Variable name "<< var_name <<" not in the bmi output vars list. \n";
      throw std::runtime_error(errMsg.str());

    }
    
  }

  // Test BMI: MODEL GRID FUNCTIONS
  std::cout<<"\n \n**************** TEST BMI GRID FUNCTIONS***********************\n";
  int grid_id[] = {0,1,2};
  int grid_size_test[] = {1,1,nz};
  int grid_rank, grid_size;
  std::string grid_type;

  for (int i=0; i< 3; i++) {
    if (VERBOSITY)  
      std::cout<<"Grid id "<< grid_id[i] <<"\n";

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_grid_rank()
    grid_rank = model.GetGridRank(grid_id[i]);
    if (grid_rank == FAILURE) return FAILURE;
    if (VERBOSITY)
      std::cout<<" rank: "<<grid_rank<<"\n";

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_grid_size
    grid_size = model.GetGridSize(grid_id[i]);
    if (grid_size == grid_size_test[i]) {
      test_status &= true;
      if (VERBOSITY)
	std::cout<<" grid size: "<<grid_size<<"\n";
    }
    else {
      test_status &= false;
      std::string passed = test_status == true ? "Yes" : "No";
      std::cout<<"Test passed: "<<passed<<"\n";
      std::stringstream errMsg;
      errMsg << "Grid size of should be "<<nz<<"\n";
      throw std::runtime_error(errMsg.str());
    }
  }
  
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /*
  // Test get_grid_type
  grid_type = model.GetGridType(grid_id);
  if (VERBOSITY)
  std::cout<<" grid type: "<<grid_type<<"\n";  
  std::cout<<"Test status: "<<test_status<<"\n";
  */

  std::cout<<GREEN<<"\n";
  std::string passed = test_status > 0 ? "Yes" : "No";
  std::cout<<"\n| *************************************** \n";
  std::cout<<"| All tests passed until this point: "<<passed<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<RESET<<"\n";

  // Test BMI: GET VALUE FUNCTIONS
  std::cout<<"\n\n************** TEST BMI GETTER SETTER FUNCTIONS********************************\n";
  
  std::cout<<"********** Input variables ***************** \n";
  // Loop through both input and output variables and call get/set_value_*()
  
  for (int i=0; i<count_in; i++) { 
    
    std::string var_name = names_in[i];
    std::cout<<"Variable name: "<< var_name <<"\n";

    // excludes the num_wetting_fronts variable which is of type int
    if (var_name != "num_wetting_fronts") {
      double *var = new double[1];
      double *dest = new double[1];
      int indices[] = {0};
      int len = 1;
      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
      // Test get_value() at each timestep
      model.GetValue(var_name, &(var[0]));
      std::cout<<" Get value: "<< var[0] <<"\n";
      
      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
      // Test get_value_at_indices()
      model.GetValueAtIndices(var_name, dest, indices, len);
      std::cout<<" Get value at indices: " << dest[0]<<"\n";
      
      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
      // Test get_value_ptr()
      double *var_ptr = new double[1];
      var_ptr= (double*) model.GetValuePtr(var_name);
      std::cout<<" Get value ptr: "<<*var_ptr<<"\n";
      
      /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
      // Test BMI set_value_at_indices()
      double dest_new[] = {-0.000472};
      
      if (var_name == "soil_storage")
	dest_new[0] = 0.8;
      
      double *dest_new_up = new double[1];
      
      //      int indices[] = {0,1,2,3}; 
      model.SetValueAtIndices(var_name, &indices[0], len, &dest_new[0]);
      
      std::cout<<" Set value at indices: "<<dest_new[0]<<"\n";
      // get_value_at_indices to see if changed
      model.GetValueAtIndices(var_name, dest_new_up,  &indices[0], len);
      std::cout<<" Get value at indices: "<<dest_new_up[0]<<"\n";
      if (dest_new[0] == dest_new_up[0])
	test_status &= true;
      else
	test_status &= false;
      
    }

  }

  passed = test_status > 0 ? "Yes" : "No";
  std::cout<<GREEN<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<"| All tests passed until this point: "<<passed<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<RESET<<"\n";
  
  std::cout<<"************* Output variables ***************** \n";

  // Update the model for Conceptual soil reservoir's unitest
  model.Update();
  
  double water_table_test_value = 1.50956; // depth from the surface in meters
  double soil_moisture_fraction_test_value = 0.174686; // soil moisture fraction in m
  double water_table_computed = 0.0;
  double soil_moisture_fraction_computed = 0.0;
    
  for (int i=0; i<count_out; i++) {
    
    std::string var_name = names_out[i];
    std::cout<<"variable name: "<< var_name <<" "<<test_status<<"\n";

    if (var_name.compare("soil_moisture_profile") == 0) {
      double *var_ptr = new double[nz];
      //var_ptr = (double*) model.GetValuePtr(var_name);
      model.GetValue(var_name, var_ptr);
      
      for (int i=0; i<nz; i++)
	std::cout<<" Get value: "<<var_ptr[i]<<"\n";

    }
    else {
      
      double val;
      
      model.GetValue(var_name, &val);
      std::cout<<" Get value: "<< val <<"\n";
      
      if (var_name.compare("soil_water_table") == 0) {
	water_table_computed = val;
	if (fabs(water_table_computed - water_table_test_value) < 0.01)
	  test_status &= true;
	else
	  test_status &= false;
      }
      else if (var_name.compare("soil_moisture_fraction") == 0) {
	soil_moisture_fraction_computed = val;
	if (fabs(soil_moisture_fraction_computed - soil_moisture_fraction_test_value) < 0.01)
	  test_status &= true;
	else
	  test_status &= false;
      }
      
    }
    
      
  }
  
  passed = test_status > 0 ? "Yes" : "No";
  std::cout<<GREEN<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<"| All BMI Tests passed: "<< RED<< passed<<RESET<<GREEN<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<RESET<<"\n";


  std::cout<<GREEN<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<"| Conceptual Soil Reservoir unittest \n";
  std::cout<<"| Unittest passed : "<< RED << passed << RESET << GREEN <<"\n";
  std::cout<<"| Water table depth (from surface) [m]: \n| Benchmark vs Computed | "
	   << water_table_test_value <<" vs "<< water_table_computed <<"\n";
  std::cout<<"| Soil moisture fraction [-]: \n| Benchmark vs Computed | "
	   << soil_moisture_fraction_test_value <<" vs "<< soil_moisture_fraction_computed <<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<RESET<<"\n";

  assert (passed == "Yes");

  //##########################################################################################################
  // Update the model for Conceptual soil reservoir's unitest

  int num_wf = 4;
  int *num_wf_ptr = &num_wf;
  
  std::string vname_smc_wf   = "soil_moisture_wetting_fronts";    // input
  std::string vname_depth_wf = "soil_depth_wetting_fronts";       // input
  std::string vname_num_wf   = "num_wetting_fronts";              // input
  std::string vname_smc_profile   = "soil_moisture_profile";      // output
  std::string vname_water_table   = "soil_water_table";           // output
  std::string vname_soil_moist_fract = "soil_moisture_fraction";  // output

  //##########################################################################################################
  // Test #1: four wetting fronts with lower two layers fully saturated; and two wetting fronts in the 2nd layer
  // benchmark values for Layered-based model
  double water_table_benchmark            = 0.44; // depth from the surface in meters
  double soil_moisture_fraction_benchmark = 0.117867; // soil moisture fraction in m
  
  {
    double soil_moisture_wetting_fronts[] = {0.25, 0.4773, 0.4773, 0.4617};
    double soil_depth_wetting_fronts[]    = {0.44, 1.12, 1.75, 2.0};
  
    // Set values
    model_layered.SetValue(vname_num_wf,num_wf_ptr);
    model_layered.SetValue(vname_smc_wf,&soil_moisture_wetting_fronts[0]);
    model_layered.SetValue(vname_depth_wf,&soil_depth_wetting_fronts[0]);
    
    model_layered.Update();
    
    water_table_computed = 0.0;
    soil_moisture_fraction_computed = 0.0;
    
    model_layered.GetValue(vname_water_table, &water_table_computed);
    model_layered.GetValue(vname_soil_moist_fract, &soil_moisture_fraction_computed);
    
    test_status = fabs(water_table_computed - water_table_benchmark) < 1.E-4 ? true : false;
    test_status &= fabs(soil_moisture_fraction_computed - soil_moisture_fraction_benchmark) < 1.E-4 ? true : false;
    
    passed = test_status > 0 ? "Yes" : "No";
    
  }
  
  std::cout<<GREEN<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<"| Layered Soil Reservoir unittest (#1) \n";
  std::cout<<"| Unittest passed : "<< RED << passed << RESET << GREEN <<"\n";
  std::cout<<"| Water table depth (from surface) [m]: \n| Benchmark vs Computed | "
	   << water_table_benchmark <<" vs "<< water_table_computed <<"\n";
  std::cout<<"| Soil moisture fraction [-]: \n| Benchmark vs Computed | "
	   << soil_moisture_fraction_benchmark <<" vs "<< soil_moisture_fraction_computed <<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<RESET<<"\n";

  assert (passed == "Yes");
 
  //##########################################################################################################
  // Test #2: four wetting fronts with only top wetting fully saturated; two wetting fronts in the top layer
  // benchmark values for Layered-based model with CONSTANT PROFILE OPTION
  water_table_benchmark            = 12.455; // depth from the surface in meters
  soil_moisture_fraction_benchmark = 0.235212; // soil moisture fraction in m
  
  {
    double soil_moisture_wetting_fronts[] = {0.4513, 0.25, 0.3173, 0.2417};
    double soil_depth_wetting_fronts[]    = {0.23, 0.4, 1.75, 2.0};
    
    // Set values
    model_layered.SetValue(vname_num_wf,num_wf_ptr);
    model_layered.SetValue(vname_smc_wf,&soil_moisture_wetting_fronts[0]);
    model_layered.SetValue(vname_depth_wf,&soil_depth_wetting_fronts[0]);
    
    model_layered.Update();
    
    water_table_computed             = 0.0;
    soil_moisture_fraction_computed = 0.0;
    
    model_layered.GetValue(vname_water_table, &water_table_computed);
    model_layered.GetValue(vname_soil_moist_fract, &soil_moisture_fraction_computed);
    
    test_status = fabs(water_table_computed - water_table_benchmark) < 1.E-4 ? true : false;
    test_status &= fabs(soil_moisture_fraction_computed - soil_moisture_fraction_benchmark) < 1.E-4 ? true : false;
    
    passed = test_status > 0 ? "Yes" : "No";
  }
  
  std::cout<<GREEN<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<"| Layered Soil Reservoir unittest (#2) \n";
  std::cout<<"| Unittest passed : "<< RED << passed << RESET << GREEN <<"\n";
  std::cout<<"| Water table depth (from surface) [m]: \n| Benchmark vs Computed | "<<water_table_benchmark<<" vs "<<water_table_computed<<"\n";
  std::cout<<"| Soil moisture fraction [-]: \n| Benchmark vs Computed | "<<soil_moisture_fraction_benchmark<<" vs "<<soil_moisture_fraction_computed<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<RESET<<"\n";
  
  assert (passed == "Yes");

  //##########################################################################################################
  // Test #3: four wetting fronts with only top wetting fully saturated; two wetting fronts in the top layer
  // benchmark values for Layered-based model with LINEAR PROFILE OPTION
  water_table_benchmark            = 12.455; // depth from the surface in meters
  soil_moisture_fraction_benchmark = 0.251803; // soil moisture fraction in m
  
  {
    model_layered.Initialize(argv[3]);
    double soil_moisture_wetting_fronts[] = {0.4513, 0.25, 0.3173, 0.2417};
    double soil_depth_wetting_fronts[]    = {0.23, 0.4, 1.75, 2.0};
    
    // Set values
    model_layered.SetValue(vname_num_wf,num_wf_ptr);
    model_layered.SetValue(vname_smc_wf,&soil_moisture_wetting_fronts[0]);
    model_layered.SetValue(vname_depth_wf,&soil_depth_wetting_fronts[0]);
    
    model_layered.Update();
    
    water_table_computed            = 0.0;
    soil_moisture_fraction_computed = 0.0;
    
    model_layered.GetValue(vname_water_table, &water_table_computed);
    model_layered.GetValue(vname_soil_moist_fract, &soil_moisture_fraction_computed);
    
    test_status = fabs(water_table_computed - water_table_benchmark) < 1.E-4 ? true : false;
    test_status &= fabs(soil_moisture_fraction_computed - soil_moisture_fraction_benchmark) < 1.E-4 ? true : false;
    
    passed = test_status > 0 ? "Yes" : "No";
  }
  
  std::cout<<GREEN<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<"| Layered Soil Reservoir unittest (#3) \n";
  std::cout<<"| Unittest passed : "<< RED << passed << RESET << GREEN <<"\n";
  std::cout<<"| Water table depth (from surface) [m]: \n| Benchmark vs Computed | "<<water_table_benchmark<<" vs "<<water_table_computed<<"\n";
  std::cout<<"| Soil moisture fraction [-]: \n| Benchmark vs Computed | "<<soil_moisture_fraction_benchmark<<" vs "<<soil_moisture_fraction_computed<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<RESET<<"\n";

  assert (passed == "Yes");

  //##########################################################################################################
  // Test #4: calibratable parameters
  {
    std::cout<<GREEN<<"\n";
    std::cout<<"| *************************************** \n";
    std::cout<<"| Calibrated parameters test \n";
    
    model.Initialize(argv[1]);

    double *smcmax_set, b_set, satpsi_set;
    double *smcmax_get, b_get, satpsi_get;

    smcmax_set = new double[1];
    smcmax_get = new double[1];
    
    // Get initial values
    model.GetValue("smcmax",&smcmax_set);
    model.GetValue("b",&b_set);
    model.GetValue("satpsi",&satpsi_set);

    std::cout<<"| Initial values | smcmax = "<< smcmax_set[0] <<", b = "<< b_set <<", satpsi = "<< satpsi_set <<"\n";
    
    // set values
    smcmax_set[0] += 0.012;
    b_set += 0.013;
    satpsi_set += 0.014;

    std::cout<<"| Setting | smcmax = "<< smcmax_set[0] <<", b = "<< b_set <<", satpsi = "<< satpsi_set <<"\n";

    model.SetValue("smcmax",&smcmax_set);
    model.SetValue("b",&b_set);
    model.SetValue("satpsi",&satpsi_set);

    model.GetValue("smcmax", &smcmax_get);
    model.GetValue("b", &b_get);
    model.GetValue("satpsi", &satpsi_get);
    
    std::cout<<"| Getting | smcmax = "<< smcmax_get[0] <<", b = "<< b_get <<", satpsi = "<< satpsi_get <<"\n";
    model.Update();
    

    test_status = fabs(smcmax_set[0] -  smcmax_get[0]) < 1.E-10 ? true : false;
    test_status &= fabs(b_set - b_get) < 1.E-10 ? true : false;
    test_status &= fabs(satpsi_set - satpsi_get) < 1.E-10 ? true : false;
    
    passed = test_status > 0 ? "Yes" : "No";

   
    std::cout<<"| Unittest passed : "<< RED << passed << RESET << GREEN <<"\n";
    std::cout<<"| *************************************** \n";
    std::cout<<RESET<<"\n";
    
    assert (passed == "Yes");
  }  

  
  return FAILURE;
}
