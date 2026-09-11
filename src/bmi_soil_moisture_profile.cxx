#ifndef BMI_SMP_C_INCLUDED
#define BMI_SMP_C_INCLUDED


#include <stdio.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "../bmi/bmi.hxx"
#include "../include/bmi_soil_moisture_profile.hxx"
#include "../include/soil_moisture_profile.hxx"
#include "Logger.hxx"

#include <boost/serialization/serialization.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

void BmiSoilMoistureProfile::
Initialize (std::string config_file)
{
#ifdef SMP_USE_EWTS    
  // Initialize the Error and Warning Trapping System
  #pragma message("SoilMoistureProfiles.bmi_soil_moisture_profile.Initialize: SMP_USE_EWTS ON")
  EwtsInit(SMP_MODULE_ID, true);
#else
  #pragma message("SoilMoistureProfiles.bmi_soil_moisture_profile.Initialize: SMP_USE_EWTS OFF")
#endif  

  LOG(LogLevel::INFO, "Initializing SMP");
  if (config_file.compare("") != 0 ) {
    this->state = new soil_moisture_profile::soil_profile_parameters;
    soil_moisture_profile::SoilMoistureProfile(config_file, state);
  }
  else {
    LOG(LogLevel::FATAL, "SMP %s config file not provided", config_file.c_str());
    throw std::runtime_error("Missing SMP Config file");
  }

  this->verbosity = this->state->verbosity;
}


void BmiSoilMoistureProfile::
Update()
{
  soil_moisture_profile::SoilMoistureProfileUpdate(state);
}


void BmiSoilMoistureProfile::
UpdateUntil(double t)
{
  soil_moisture_profile::SoilMoistureProfileUpdate(state);
}


void BmiSoilMoistureProfile::
Finalize()
{
  if (this->state) {
    delete state;
    this->state = NULL;
  }
}

void BmiSoilMoistureProfile::
PrintSoilMoistureProfile()
{
  soil_moisture_profile::PrintSoilMoistureProfile(this->state);
}

int BmiSoilMoistureProfile::
GetVarGrid(std::string name)
{
  if (
    name.compare("soil_storage_model") == 0
    || name.compare("num_wetting_fronts") == 0
    || name.compare("serialization_free") == 0
  ) // int
    return 0;
  else if (
    name.compare("soil_storage") == 0
    || name.compare("soil_storage_change") == 0
	  || name.compare("soil_water_table") == 0
    || name.compare("soil_moisture_fraction") == 0
    || name.compare("Qb_topmodel") == 0
    || name.compare("Qv_topmodel") == 0
    || name.compare("global_deficit") == 0
    || name.compare("b") == 0
    || name.compare("satpsi") == 0
    || name.compare("reset_time") == 0
  ) // double
    return 1;
  else if (
    name.compare("soil_moisture_profile") == 0
  ) // array of doubles (conceptual model)
    return 2;
  else if (
    name.compare("soil_moisture_wetting_fronts") == 0
    || name.compare("soil_depth_wetting_fronts") == 0
  ) // array of doubles (layered model)
    return 3;
  else if (
    name.compare("smcmax") == 0
  ) // fixed number of layers for calibratable params
    return 4;
  else if (
    name.compare("serialization_state") == 0
  ) // char* beginning of serialized data
    return 5;
  else if (
    name.compare("serialization_create") == 0
    || name.compare("serialization_size") == 0
  ) // uint64_t
    return 6;
  else
    return -1;
}


std::string BmiSoilMoistureProfile::
GetVarType(std::string name)
{
  int var_grid = GetVarGrid(name);

  if (var_grid == 0)
    return "int";
  else if (var_grid == 1 || var_grid == 2 || var_grid == 3 || var_grid == 4)
    return "double";
  else if (var_grid == 5)
    return "char";
  else if (var_grid == 6)
    return "uint64_t";
  else
    return "";
}


int BmiSoilMoistureProfile::
GetVarItemsize(std::string name)
{
  std::string var_type = GetVarType(name);

  if (var_type.compare("int") == 0)
    return sizeof(int);
  else if (var_type.compare("double") == 0)
    return sizeof(double);
  else if (var_type.compare("char") == 0)
    return sizeof(char);
  else if (var_type.compare("uint64_t") == 0)
    return sizeof(uint64_t);
  else
    return 0;
}


std::string BmiSoilMoistureProfile::
GetVarUnits(std::string name)
{
  if (name.compare("soil_storage") == 0 || name.compare("soil_storage_change") == 0 ||
      name.compare("soil_water_table") == 0)
    return "m";
  else if (name.compare("soil_moisture_profile") == 0 || name.compare("soil_moisture_wetting_fronts") == 0 ||
	   name.compare("soil_moisture_fraction") == 0)
    return "1";  // dimensionless (UDUNITS)
  else if (name.compare("Qb_topmodel") == 0 || name.compare("Qv_topmodel") == 0)
    return "m h^-1";
  else if (name.compare("global_deficit") == 0)
    return "m";
  else if (name.compare("soil_depth_wetting_fronts") == 0)
    return "m";
  else
    return "none";
}


int BmiSoilMoistureProfile::
GetVarNbytes(std::string name)
{
  int itemsize;
  int gridsize;

  itemsize = this->GetVarItemsize(name);
  gridsize = this->GetGridSize(this->GetVarGrid(name));
  return itemsize * gridsize;
}


std::string BmiSoilMoistureProfile::
GetVarLocation(std::string name)
{
  int var_grid = GetVarGrid(name);

  if (var_grid <= 3 && name != "serialization_free")
    return "node";
  else
    return "none";
}


void BmiSoilMoistureProfile::
GetGridShape(const int grid, int *shape)
{
  if (grid == 2) {
    shape[0] = this->state->shape[0];
  }
  else if (grid == 3) {
    shape[0] = this->state->shape[1];
  }
  else if (grid == 4) {
    shape[0] = this->state->num_layers;
  }
}


void BmiSoilMoistureProfile::
GetGridSpacing (const int grid, double * spacing)
{
  if (grid == 0) {
    spacing[0] = this->state->spacing[0];
  }
}


void BmiSoilMoistureProfile::
GetGridOrigin (const int grid, double *origin)
{
  if (grid == 0) {
    origin[0] = this->state->origin[0];
  }
}


int BmiSoilMoistureProfile::
GetGridRank(const int grid)
{
  if (grid <= 3)
    return 1;
  else
    return -1;
}


int BmiSoilMoistureProfile::
GetGridSize(const int grid)
{
  if (grid == 0 || grid == 1 || grid == 6)
    return 1;
  else if (grid == 2)
    return this->state->shape[0];
  else if (grid == 3)
    return this->state->shape[1];
  else if (grid == 4)
    return this->state->num_layers;
  else if (grid == 5)
    return this->m_serialized_length; // size of currently saved state
  else
    return -1;
}


std::string BmiSoilMoistureProfile::
GetGridType(const int grid)
{
  if (grid == 0)
    return "uniform_rectilinear";
  else
    return "";
}


void BmiSoilMoistureProfile::
GetGridX(const int grid, double *x)
{
  LOG(LogLevel::SEVERE, "GetGridX Not Implemented");
  throw coupler::NotImplemented();
}


void BmiSoilMoistureProfile::
GetGridY(const int grid, double *y)
{
  LOG(LogLevel::SEVERE, "GetGridY Not Implemented");
  throw coupler::NotImplemented();
}


void BmiSoilMoistureProfile::
GetGridZ(const int grid, double *z)
{
  LOG(LogLevel::SEVERE, "GetGridZ Not Implemented");
  throw coupler::NotImplemented();
}


int BmiSoilMoistureProfile::
GetGridNodeCount(const int grid)
{
  LOG(LogLevel::SEVERE, "GetGridNodeCount Not Implemented");
  throw coupler::NotImplemented();
  /*
  if (grid == 0)
    return this->state->shape[0];
  else
    return -1;
  */
}


int BmiSoilMoistureProfile::
GetGridEdgeCount(const int grid)
{
  LOG(LogLevel::SEVERE, "GetGridEdgeCount Not Implemented");
  throw coupler::NotImplemented();
}


int BmiSoilMoistureProfile::
GetGridFaceCount(const int grid)
{
  LOG(LogLevel::SEVERE, "GetGridFaceCount Not Implemented");
  throw coupler::NotImplemented();
}


void BmiSoilMoistureProfile::
GetGridEdgeNodes(const int grid, int *edge_nodes)
{
  LOG(LogLevel::SEVERE, "GetGridEdgeNodes Not Implemented");
  throw coupler::NotImplemented();
}


void BmiSoilMoistureProfile::
GetGridFaceEdges(const int grid, int *face_edges)
{
  LOG(LogLevel::SEVERE, "GetGridFaceEdges Not Implemented");
  throw coupler::NotImplemented();
}


void BmiSoilMoistureProfile::
GetGridFaceNodes(const int grid, int *face_nodes)
{
  LOG(LogLevel::SEVERE, "GetGridFaceNodes Not Implemented");
  throw coupler::NotImplemented();
}


void BmiSoilMoistureProfile::
GetGridNodesPerFace(const int grid, int *nodes_per_face)
{
  LOG(LogLevel::SEVERE, "GetGridNodesPerFace Not Implemented");
  throw coupler::NotImplemented();
}


void BmiSoilMoistureProfile::
GetValue (std::string name, void *dest)
{
  void * src = NULL;
  int nbytes = 0;

  src = this->GetValuePtr(name);
  nbytes = this->GetVarNbytes(name);
  memcpy (dest, src, nbytes);
}


void *BmiSoilMoistureProfile::
GetValuePtr (std::string name)
{
  if (name.compare("soil_storage") == 0)
    return (void*)(&this->state->soil_storage);
  else if (name.compare("soil_storage_change") == 0)
    return (void*)(&this->state->soil_storage_change_per_timestep);
  else  if (name.compare("soil_water_table") == 0)
    return (void*)(&this->state->water_table_depth);
  else  if (name.compare("soil_moisture_fraction") == 0)
    return (void*)(&this->state->soil_moisture_fraction);
  else if (name.compare("soil_moisture_profile") == 0)
    return (void*)this->state->soil_moisture_profile;
  else if (name.compare("soil_moisture_wetting_fronts") == 0)
    return this->state->soil_moisture_wetting_fronts.data();
  else if (name.compare("soil_depth_wetting_fronts") == 0)
    return this->state->soil_depth_wetting_fronts.data();
  else if (name.compare("soil_storage_model") == 0)
    return (void*)(&this->state->soil_storage_model);
  else if (name.compare("num_wetting_fronts") == 0)
    return (void*)(&this->state->num_wetting_fronts);
  else if (name.compare("Qb_topmodel") == 0)
    return (void*)(&this->state->Qb_topmodel);
  else if (name.compare("Qv_topmodel") == 0)
    return (void*)(&this->state->Qv_topmodel);
  else if (name.compare("global_deficit") == 0)
    return (void*)(&this->state->global_deficit);
  else if (name.compare("smcmax") == 0)
    return (void*)this->state->smcmax;
  else if (name.compare("b") == 0)
    return (void*)(&this->state->b);
  else if (name.compare("satpsi") == 0)
    return (void*)(&this->state->satpsi);
  else if (name.compare("serialization_state") == 0)
    return (void*)(this->m_serialized.data());
  else if (name.compare("serialization_size") == 0) {
    return (void*)(&this->m_serialized_length);
  } else {
    std::stringstream errMsg;
    errMsg << "variable "<< name << " does not exist";
    LOG(LogLevel::FATAL, errMsg.str());
    throw std::runtime_error(errMsg.str());
    return NULL;
  }
}


void BmiSoilMoistureProfile::
GetValueAtIndices (std::string name, void *dest, int *inds, int len)
{
  void * src = NULL;

  src = this->GetValuePtr(name);

  if (src) {
    int i;
    int itemsize = 0;
    int offset;
    char *ptr;

    itemsize = this->GetVarItemsize(name);
    for (i=0, ptr=(char *)dest; i<len; i++, ptr+=itemsize) {
      offset = inds[i] * itemsize;
      memcpy(ptr, (char *)src + offset, itemsize);
    }
  }
}

void BmiSoilMoistureProfile::
ResetSize (std::string name)
{
// reset the size of wetting fronts array to the number of wetting fronts at the timestep
  if (name.compare("soil_moisture_wetting_fronts") == 0) {
    if (this->state->num_wetting_fronts <= 0) {
      std::string error_msg = "The number of wetting fronts must be greater than zero. The current number of wetting fronts is " + std::to_string(this->state->num_wetting_fronts);
      LOG(LogLevel::FATAL, error_msg);
      throw std::out_of_range(error_msg);
    }
    state->soil_moisture_wetting_fronts.resize(this->state->num_wetting_fronts);
  }
  else if (name.compare("soil_depth_wetting_fronts") == 0) {
    if (this->state->num_wetting_fronts <= 0) {
      std::string error_msg = "The number of wetting fronts must be greater than zero. The current number of wetting fronts is " + std::to_string(this->state->num_wetting_fronts);
      LOG(LogLevel::FATAL, error_msg);
      throw std::out_of_range(error_msg);
    }
    state->soil_depth_wetting_fronts.resize(this->state->num_wetting_fronts);
  }
}

void BmiSoilMoistureProfile::
SetValue (std::string name, void *src)
{
  // special cases for state serialization
  if (name.compare("serialization_state") == 0) {
    this->load_serialized((char*)src);
    return;
  } else if (name.compare("serialization_create") == 0) {
    this->new_serialized();
    return;
  } else if (name.compare("serialization_free") == 0) {
    this->free_serialized();
    return;
  } else if (name == "reset_time") {
    // This BMI does not increment its time, so no action needs to be done.
    return;
  }
  void * dest = NULL;
  ResetSize(name);
  
  dest = this->GetValuePtr(name);
  
  if (dest) {
    int nbytes = 0;
    nbytes = this->GetVarNbytes(name);
    memcpy(dest, src, nbytes);
    
    if (name.compare("num_wetting_fronts") == 0)
      this->state->shape[1] = this->state->num_wetting_fronts;
    
  }
  
}


void BmiSoilMoistureProfile::
SetValueAtIndices (std::string name, int * inds, int len, void *src)
{
  void * dest = NULL;

  ResetSize(name);
  
  dest = this->GetValuePtr(name);
  
  if (dest) {
    int i;
    int itemsize = 0;
    int offset;
    char *ptr;

    itemsize = this->GetVarItemsize(name);
    
    for (i=0, ptr=(char *)src; i<len; i++, ptr+=itemsize) {
      offset = inds[i] * itemsize;
      memcpy((char *)dest + offset, ptr, itemsize);

      if (name.compare("num_wetting_fronts") == 0)
	this->state->shape[1] = this->state->num_wetting_fronts;
    }
  }
}


std::string BmiSoilMoistureProfile::
GetComponentName()
{
  return "SoilMoistureProfiles BMI";
}


int BmiSoilMoistureProfile::
GetInputItemCount()
{
  return this->input_var_name_count;
}


int BmiSoilMoistureProfile::
GetOutputItemCount()
{
  return this->output_var_name_count;
}


std::vector<std::string> BmiSoilMoistureProfile::
GetInputVarNames()
{
  std::vector<std::string> names;
  
  for (int i=0; i<this->input_var_name_count; i++)
    names.push_back(this->input_var_names[i]);
  
  return names;
}


std::vector<std::string> BmiSoilMoistureProfile::
GetOutputVarNames()
{
  std::vector<std::string> names;

  for (int i=0; i<this->output_var_name_count; i++)
    names.push_back(this->output_var_names[i]);

  return names;
}


double BmiSoilMoistureProfile::
GetStartTime () {
  return 0.0;
}


double BmiSoilMoistureProfile::
GetEndTime () {
  return 0.0;
}


double BmiSoilMoistureProfile::
GetCurrentTime () {
  return 0.0;
}


std::string BmiSoilMoistureProfile::
GetTimeUnits() {
  return "s";
}


double BmiSoilMoistureProfile::
GetTimeStep () {
  return 0;
}


template<class Archive>
void BmiSoilMoistureProfile::
serialize(Archive& ar, const unsigned int version) {
  soil_moisture_profile::soil_profile_parameters* state = this->state;
  // size of array pointers assigned in initialization
  int size = state->ncells;

  // all three models (Conceptual, Layered, and Topmodel) create these three states
  ar & state->init_profile;
  ar & boost::serialization::make_array(state->soil_moisture_profile, size);
  ar & state->water_table_depth;

  // state regardless of model
  ar & state->soil_storage;
  ar & state->soil_moisture_fraction;
}


void BmiSoilMoistureProfile::
new_serialized() {
  // resize with space for size as a header
  this->m_serialized.resize(sizeof(uint64_t));
  boost::archive::binary_oarchive archive(this->m_serialized);
  try {
    archive << (*this);
    this->m_serialized_length = this->m_serialized.size();
    // store size of serialized data as header
    uint64_t serialized_size = this->m_serialized_length - sizeof(uint64_t);
    memcpy(this->m_serialized.data(), &serialized_size, sizeof(uint64_t));
  } catch (const std::exception &e) {
    LOG(LogLevel::WARNING, "Serializing SMP encounterd an error: %s", e.what());
    LOG(LogLevel::WARNING, "Set m_serialized_length = 0");
    this->m_serialized_length = 0;
    throw;
  }
}


void BmiSoilMoistureProfile::
load_serialized(char* data) {
  // get size from header of data
  uint64_t size;
  memcpy(&size, data, sizeof(uint64_t));
  // create stream from after header
  membuf stream(data + sizeof(uint64_t), size);
  boost::archive::binary_iarchive archive(stream);
  try {
    archive >> (*this);
  } catch (const std::exception &e) {
    LOG(LogLevel::SEVERE, "Deserializing SMP encounterd an error: %s", e.what());
    throw;
  }
  this->free_serialized();
}


void BmiSoilMoistureProfile::
free_serialized() {
  this->m_serialized.clear();
  this->m_serialized.shrink_to_fit();
  this->m_serialized_length = 0;
}


#endif
