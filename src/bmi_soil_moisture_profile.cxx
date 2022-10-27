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


void BmiSoilMoistureProfile::
Initialize (std::string config_file)
{
  if (config_file.compare("") != 0 ) {
    this->model = new soil_moisture_profile::soil_profile_parameters;
    soil_moisture_profile::SoilMoistureProfile(config_file, model);
  }
}


void BmiSoilMoistureProfile::
Update()
{
  SoilMoistureProfileUpdate(model);
}


void BmiSoilMoistureProfile::
UpdateUntil(double t)
{
  SoilMoistureProfileUpdate(model);
}


void BmiSoilMoistureProfile::
Finalize()
{
  if (this->model)
    delete model;
}


int BmiSoilMoistureProfile::
GetVarGrid(std::string name)
{
  if (name.compare("soil_storage_model") == 0 || name.compare("num_cells_layered") == 0)   // int
    return 0;
  else if (name.compare("soil_storage") == 0 || name.compare("soil_storage_change") == 0 || name.compare("soil_water_table") == 0) // double
    return 1; 
  else if (name.compare("soil_moisture_profile") == 0) // array of doubles (conceptual model)
    return 2;
  else if (name.compare("soil_moisture_layered") == 0 || name.compare("soil_depths_layered") == 0) // array of doubles (layered model)
    return 3; 
  else
    return -1;
}


std::string BmiSoilMoistureProfile::
GetVarType(std::string name)
{
  if (name.compare("soil_storage_model") == 0 || name.compare("num_cells_layered") == 0)
    return "int";
  else if (name.compare("soil_storage") == 0 || name.compare("soil_storage_change") == 0 || name.compare("soil_water_table") == 0)
    return "double";
  else if (name.compare("soil_moisture_profile") == 0 || name.compare("soil_moisture_layered") == 0 || name.compare("soil_depths_layered") == 0)
    return "double";
  else
    return "";
}


int BmiSoilMoistureProfile::
GetVarItemsize(std::string name)
{
  if (name.compare("soil_storage_model") == 0 || name.compare("num_cells_layered") == 0)
    return sizeof(int);
  else if (name.compare("soil_storage") == 0 || name.compare("soil_storage_change") == 0 || name.compare("soil_water_table") == 0)
    return sizeof(double);
  else if (name.compare("soil_moisture_profile") == 0 || name.compare("soil_moisture_layered") == 0 || name.compare("soil_depths_layered") == 0)
    return sizeof(double);
  else
    return 0;
}


std::string BmiSoilMoistureProfile::
GetVarUnits(std::string name)
{
  if (name.compare("soil_storage") == 0 || name.compare("soil_storage_change") == 0  || name.compare("soil_water_table") == 0)
    return "m";
  else if (name.compare("soil_moisture_profile") == 0 || name.compare("soil_moisture_layered") == 0)
    return "none";
  else if (name.compare("soil_depths_layered") == 0)
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
  if (name.compare("soil_storage") == 0 || name.compare("soil_storage_change") == 0 || name.compare("soil_water_table") == 0)
    return "node";
  else if (name.compare("soil_moisture_profile") == 0 || name.compare("soil_moisture_layered") == 0 || name.compare("soil_depths_layered") == 0 || name.compare("num_cells_layered") == 0)
    return "node";
  else
    return "none";
}


void BmiSoilMoistureProfile::
GetGridShape(const int grid, int *shape)
{
  if (grid == 2) {
    shape[0] = this->model->shape[0];
  }
  else if (grid == 3) {
    shape[0] = this->model->shape[1];
  }
}


void BmiSoilMoistureProfile::
GetGridSpacing (const int grid, double * spacing)
{
  if (grid == 0) {
    spacing[0] = this->model->spacing[0];
  }
}


void BmiSoilMoistureProfile::
GetGridOrigin (const int grid, double *origin)
{
  if (grid == 0) {
    origin[0] = this->model->origin[0];
  }
}


int BmiSoilMoistureProfile::
GetGridRank(const int grid)
{
  if (grid == 0 || grid == 1 || grid == 2)
    return 1;
  else
    return -1;
}


int BmiSoilMoistureProfile::
GetGridSize(const int grid)
{
  if (grid == 0 || grid == 1)
    return 1;
  else if (grid == 2)
    return this->model->shape[0];
  else if (grid == 3)
    return this->model->shape[1];
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
  /*
  if (grid == 0)
    return this->model->shape[0];
  else
    return -1;
  */
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
    return (void*)(&this->model->soil_storage);
  else if (name.compare("soil_storage_change") == 0)
    return (void*)(&this->model->soil_storage_change_per_timestep);
  else  if (name.compare("soil_water_table") == 0)
    return (void*)(&this->model->water_table_depth);
  else if (name.compare("soil_moisture_profile") == 0)
    return (void*)this->model->soil_moisture_profile;
  else if (name.compare("soil_moisture_layered") == 0)
    return (void*)this->model->soil_moisture_layered;
  else if (name.compare("soil_depths_layered") == 0)
    return (void*)this->model->soil_depths_layered;
  else if (name.compare("soil_storage_model") == 0)
    return (void*)(&this->model->soil_storage_model);
  else if (name.compare("num_cells_layered") == 0)
    return (void*)(&this->model->ncells_layered);
  else {
    std::stringstream errMsg;
    errMsg << "variable "<< name << " does not exist";
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
SetValue (std::string name, void *src)
{
  void * dest = NULL;
  dest = this->GetValuePtr(name);

  if (dest) {
    int nbytes = 0;
    nbytes = this->GetVarNbytes(name);
    memcpy(dest, src, nbytes);
  }

}


void BmiSoilMoistureProfile::
SetValueAtIndices (std::string name, int * inds, int len, void *src)
{
  void * dest = NULL;

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
    }
  }
}


std::string BmiSoilMoistureProfile::
GetComponentName()
{
  return "SoilMoistureProfile BMI";
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

#endif
