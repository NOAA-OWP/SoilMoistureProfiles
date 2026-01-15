#ifndef BMI_SMP_C_INCLUDED
#define BMI_SMP_C_INCLUDED


#include <string>
#include <cstring>
#include <vector>
#include "../include/bmi_soil_moisture_profile.hxx"
#include "../include/soil_moisture_profile.hxx"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

void BmiSoilMoistureProfile::
Initialize (const std::string config_file)
{
  if (config_file.empty()) {
    this->state = new soil_moisture_profile::soil_profile_parameters;
    soil_moisture_profile::SoilMoistureProfile(config_file, state);
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
    this->state = nullptr;
  }
}

void BmiSoilMoistureProfile::
PrintSoilMoistureProfile()
{
  soil_moisture_profile::PrintSoilMoistureProfile(this->state);
}

int BmiSoilMoistureProfile::
GetVarGrid(const std::string name) {
  if (
    name == "soil_storage_model"
    || name == "num_wetting_fronts"
    || name == "serialization_free"
  ) {
    return 0;
  }

  if (
    name == "soil_storage"
    || name == "soil_storage_change"
	  || name == "soil_water_table"
    || name == "soil_moisture_fraction"
    || name == "Qb_topmodel"
    || name == "Qv_topmodel"
    || name == "global_deficit"
    || name == "b"
    || name == "satpsi"
  ) {
    return 1;
  }

  if (name == "soil_moisture_profile") {
    // array of doubles (conceptual model)
    return 2;
  }
  if (name == "soil_moisture_wetting_fronts" || name == "soil_depth_wetting_fronts") {
    // array of doubles (layered model)
    return 3;
  }
  if (name == "smcmax") {
    // fixed number of layers for calibratable params
    return 4;
  }
  if (name == "serialization_state") // char* beginning of serialized data
    return 5;
  if (name == "serialization_create" || name == "serialization_size") {
    // uint64_t
    return 6;
  }

  return -1;
}


std::string BmiSoilMoistureProfile::
GetVarType(const std::string name) {
  const int var_grid = GetVarGrid(name);

  if (var_grid == 0)
    return "int";
  if (var_grid == 1 || var_grid == 2 || var_grid == 3 || var_grid == 4)
    return "double";
  if (var_grid == 5)
    return "char";
  if (var_grid == 6)
    return "uint64_t";
  return "";
}


int BmiSoilMoistureProfile::
GetVarItemsize(const std::string name)
{
  const std::string var_type = GetVarType(name);

  if (var_type == "int")
    return sizeof(int);
  if (var_type == "double")
    return sizeof(double);
  if (var_type == "char")
    return sizeof(char);
  if (var_type == "uint64_t")
    return sizeof(uint64_t);

  return 0;
}


std::string BmiSoilMoistureProfile::
GetVarUnits(std::string name)
{
  if (name == "soil_storage" || name == "soil_storage_change" ||
      name == "soil_water_table")
    return "m";
  if (name == "soil_moisture_profile" || name == "soil_moisture_wetting_fronts" ||
	   name == "soil_moisture_fraction")
    return "1";  // dimensionless (UDUNITS)
  if (name == "Qb_topmodel" || name == "Qv_topmodel")
    return "m h^-1";
  if (name == "global_deficit")
    return "m";
  if (name == "soil_depth_wetting_fronts")
    return "m";

  return "none";
}


int BmiSoilMoistureProfile::
GetVarNbytes(const std::string name)
{
  const int item_size = this->GetVarItemsize(name);;
  const int grid_size = this->GetGridSize(this->GetVarGrid(name));
  return item_size * grid_size;
}


std::string BmiSoilMoistureProfile::
GetVarLocation(const std::string name)
{
  const int var_grid = GetVarGrid(name);

  if (var_grid <= 3 && name != "serialization_free")
    return "node";

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

  return -1;
}


int BmiSoilMoistureProfile::
GetGridSize(const int grid)
{
  if (grid == 0 || grid == 1 || grid == 6)
    return 1;
  if (grid == 2)
    return this->state->shape[0];
  if (grid == 3)
    return this->state->shape[1];
  if (grid == 4)
    return this->state->num_layers;
  if (grid == 5) {
    // WARNING: The narrowing of a uint64_t to an int is implementation-defined and not guaranteed
    return this->m_serialized_length; // size of currently saved state
  }
  return -1;
}


std::string BmiSoilMoistureProfile::
GetGridType(const int grid)
{
  if (grid == 0)
    return "uniform_rectilinear";
  return "";
}


void BmiSoilMoistureProfile::
GetGridX(const int grid, double *x)
{
  throw coupler::NotImplemented();
}


void BmiSoilMoistureProfile::
GetGridY(const int grid, double *y)
{
  throw coupler::NotImplemented();
}


void BmiSoilMoistureProfile::
GetGridZ(const int grid, double *z)
{
  throw coupler::NotImplemented();
}


int BmiSoilMoistureProfile::
GetGridNodeCount(const int grid)
{
  throw coupler::NotImplemented();
}


int BmiSoilMoistureProfile::
GetGridEdgeCount(const int grid)
{
  throw coupler::NotImplemented();
}


int BmiSoilMoistureProfile::
GetGridFaceCount(const int grid)
{
  throw coupler::NotImplemented();
}


void BmiSoilMoistureProfile::
GetGridEdgeNodes(const int grid, int *edge_nodes)
{
  throw coupler::NotImplemented();
}


void BmiSoilMoistureProfile::
GetGridFaceEdges(const int grid, int *face_edges)
{
  throw coupler::NotImplemented();
}


void BmiSoilMoistureProfile::
GetGridFaceNodes(const int grid, int *face_nodes)
{
  throw coupler::NotImplemented();
}


void BmiSoilMoistureProfile::
GetGridNodesPerFace(const int grid, int *nodes_per_face)
{
  throw coupler::NotImplemented();
}


void BmiSoilMoistureProfile::
GetValue (const std::string name, void *dest)
{
  const void *src = this->GetValuePtr(name);
  int nbytes = this->GetVarNbytes(name);
  memcpy(dest, src, nbytes);
}


void *BmiSoilMoistureProfile::
GetValuePtr (const std::string name)
{
  if (name == "soil_storage")
    return &this->state->soil_storage;
  if (name == "soil_storage_change")
    return &this->state->soil_storage_change_per_timestep;
  if (name == "soil_water_table")
    return &this->state->water_table_depth;
  if (name == "soil_moisture_fraction")
    return &this->state->soil_moisture_fraction;
  if (name == "soil_moisture_profile")
    return this->state->soil_moisture_profile;
  if (name == "soil_moisture_wetting_fronts")
    return this->state->soil_moisture_wetting_fronts;
  if (name == "soil_depth_wetting_fronts")
    return this->state->soil_depth_wetting_fronts;
  if (name == "soil_storage_model")
    return &this->state->soil_storage_model;
  if (name == "num_wetting_fronts")
    return &this->state->num_wetting_fronts;
  if (name == "Qb_topmodel")
    return &this->state->Qb_topmodel;
  if (name == "Qv_topmodel")
    return &this->state->Qv_topmodel;
  if (name == "global_deficit")
    return &this->state->global_deficit;
  if (name == "smcmax")
    return this->state->smcmax;
  if (name == "b")
    return &this->state->b;
  if (name == "satpsi")
    return &this->state->satpsi;
  if (name == "serialization_state")
    return (void*)(this->m_serialized.data());
  if (name == "serialization_size") {
    return &this->m_serialized_length;
  }

  std::stringstream errMsg;
  errMsg << "variable "<< name << " does not exist";
  throw std::runtime_error(errMsg.str());
}


void BmiSoilMoistureProfile::
GetValueAtIndices (std::string name, void *dest, int *inds, int len)
{
  void * src = this->GetValuePtr(name);

  if (src) {
    int i;
    int item_size = 0;
    int offset;
    char *ptr;

    item_size = this->GetVarItemsize(name);
    for (i=0, ptr=(char *)dest; i<len; i++, ptr+=item_size) {
      offset = inds[i] * item_size;
      memcpy(ptr, (char *)src + offset, item_size);
    }
  }
}

void BmiSoilMoistureProfile::
ResetSize (const std::string name)
{
// reset the size of wetting fronts array to the number of wetting fronts at the timestep
  if (name == "soil_moisture_wetting_fronts") {
    assert (this->state->num_wetting_fronts > 0);
    state->soil_moisture_wetting_fronts = new double[this->state->num_wetting_fronts]();
  }
  else if (name == "soil_depth_wetting_fronts") {
    assert (this->state->num_wetting_fronts > 0);
    state->soil_depth_wetting_fronts = new double[this->state->num_wetting_fronts]();
  }
}

void BmiSoilMoistureProfile::
SetValue (const std::string name, void *src)
{
  // special cases for state serialization
  if (name == "serialization_state") {
    this->load_serialized(static_cast<char*>(src));
    return;
  }

  if (name == "serialization_create") {
    this->new_serialized();
    return;
  }

  if (name == "serialization_free") {
    this->free_serialized();
    return;
  }

  ResetSize(name);

  void *dest = this->GetValuePtr(name);

  if (dest) {
    int nbytes = 0;
    nbytes = this->GetVarNbytes(name);
    memcpy(dest, src, nbytes);

    if (name == "num_wetting_fronts")
      this->state->shape[1] = this->state->num_wetting_fronts;

  }

}


void BmiSoilMoistureProfile::
SetValueAtIndices (std::string name, int * inds, int len, void *src)
{
  ResetSize(name);

  void * dest = this->GetValuePtr(name);

  if (dest) {
    int i;
    const int item_size = this->GetVarItemsize(name);
    char *ptr;

    for (i=0, ptr=static_cast<char*>(src); i<len; i++, ptr += item_size) {
      const int offset = inds[i] * item_size;
      memcpy(static_cast<char*>(dest) + offset, ptr, item_size);

      if (name == "num_wetting_fronts")
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
  return input_var_name_count;
}


int BmiSoilMoistureProfile::
GetOutputItemCount()
{
  return output_var_name_count;
}


std::vector<std::string> BmiSoilMoistureProfile::
GetInputVarNames()
{
  std::vector<std::string> names;

  for (const std::string& input_var_name : this->input_var_names)
    names.push_back(input_var_name);

  return names;
}


std::vector<std::string> BmiSoilMoistureProfile::
GetOutputVarNames()
{
  std::vector<std::string> names;

  for (const std::string& output_var_name : this->output_var_names)
    names.push_back(output_var_name);

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
  std::cerr <<
    "WARNING: Serialization has not been properly implemented in BmiSoilMoistureProfile or its derivatives"
  << std::endl;
  soil_moisture_profile::soil_profile_parameters* state = this->state;
  // size of array pointers assigned in initialization
  const int size = state->ncells;

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
  std::cerr << "WARNING: Serialization has not been properly implemented in BmiSoilMoistureProfile or its derivatives" << std::endl;
  this->m_serialized.clear();
  boost::archive::binary_oarchive archive(this->m_serialized);
  try {
    archive << (*this);
    this->m_serialized_length = this->m_serialized.size();
  } catch (const std::exception &e) {
    // stringstream ss;
    // ss << "Serializing SMP encounterd an error: " << e.what();
    // Logger::Log(ss.str(), LogLevel::SEVERE);
    this->m_serialized_length = 0;
    throw;
  }
}


void BmiSoilMoistureProfile::
load_serialized(const char* data) {
  std::cerr << "WARNING: Serialization has not been properly implemented in BmiSoilMoistureProfile or its derivatives" << std::endl;
  std::stringstream stream(data);
  boost::archive::binary_iarchive archive(stream);
  try {
    archive >> (*this);
  } catch (const std::exception &e) {
    // stringstream ss;
    // ss << "Deserializing SMP encounterd an error: " << e.what();
    // Logger::Log(ss.str(), LogLevel::SEVERE);
    throw;
  }
  this->free_serialized();
}


void BmiSoilMoistureProfile::
free_serialized() {
  std::cerr << "WARNING: Serialization has not been properly implemented in BmiSoilMoistureProfile or its derivatives" << std::endl;
  this->m_serialized.clear();
  this->m_serialized.shrink_to_fit();
  this->m_serialized_length = 0;
}


#endif
