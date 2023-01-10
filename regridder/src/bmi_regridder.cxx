#ifndef BMI_MAPPER_C_INCLUDED
#define BMI_MAPPER_C_INCLUDED


#include <stdio.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include "../bmi/bmi.hxx"
#include "../include/bmi_regridder.hxx"
#include "../include/soil_moisture_regridding.hxx"

void BmiMapper::
Initialize (std::string config_file)
{
  if (config_file.compare("") != 0 )
    this->_model = new smc_mapping::SoilMoistureMapping(config_file);
}


void BmiMapper::
Update()
{
  this->_model->SoilMoistureFromBasinToGrid();
}


void BmiMapper::
UpdateUntil(double t)
{
  this->_model->SoilMoistureFromBasinToGrid();
}


void BmiMapper::
Finalize()
{
  this->_model->~SoilMoistureMapping();
}


int BmiMapper::
GetVarGrid(std::string name)
{
  if (name.compare("area_fraction") == 0 || name.compare("grid_gid") == 0 || name.compare("grid_gid_unique") == 0 || name.compare("grid_soil_moisture") == 0)
    return 0;
  else if (name.compare("TWI") == 0 || name.compare("global_deficit") == 0 || name.compare("porosity") == 0 || name.compare("depth") == 0)
    return 1;
  else 
    return -1;
}


std::string BmiMapper::
GetVarType(std::string name)
{
  if (name.compare("area_fraction") == 0 || name.compare("grid_soil_moisture") == 0)
    return "double";
  else if (name.compare("TWI") == 0 || name.compare("global_deficit") == 0 || name.compare("porosity") == 0 || name.compare("depth") == 0)
    return "double";
  else if (name.compare("grid_gid") == 0 || name.compare("grid_gid_unique") == 0)
    return "int";
  else
    return "";
}


int BmiMapper::
GetVarItemsize(std::string name)
{
  if (name.compare("area_fraction") == 0 || name.compare("grid_soil_moisture") == 0)
    return sizeof("double");
  if (name.compare("TWI") == 0 || name.compare("global_deficit") == 0 || name.compare("porosity") == 0 || name.compare("depth") == 0 )
    return sizeof(double);
  else if (name.compare("grid_gid") == 0 || name.compare("grid_gid_unique") == 0)
    return sizeof(int);
  else
    return 0;
}


std::string BmiMapper::
GetVarUnits(std::string name)
{
  if (name.compare("TWI") == 0)
    return "ln(m)";
  else if (name.compare("global_deficit") == 0)
    return "m";
  else if (name.compare("depth") == 0)
    return "m";
  else if (name.compare("area_fraction") == 0)
    return "m^-2"; //check the notation - ajk
  else
    return "none";
}


int BmiMapper::
GetVarNbytes(std::string name)
{
  int itemsize;
  int gridsize;

  itemsize = this->GetVarItemsize(name);
  gridsize = this->GetGridSize(this->GetVarGrid(name));
  return itemsize * gridsize;
}


std::string BmiMapper::
GetVarLocation(std::string name)
{
  if (name.compare("TWI") == 0 || name.compare("global_deficit") == 0 || name.compare("porosity") == 0 || name.compare("depth") == 0 || name.compare("area_fraction") == 0 || name.compare("grid_gid") == 0 || name.compare("grid_gid_unique") == 0 || name.compare("grid_soil_moisture") == 0)
    return "node";
  else
    return "";
}


void BmiMapper::
GetGridShape(const int grid, int *shape)
{
  //NWM grid
  if (grid == 0) {
    shape[0] = this->_model->shape[0];
  }
  // basins, not sure if we will need this???
  if (grid == 1) {
    shape[0] = this->_model->shape_basin[0];
  }
}


void BmiMapper::
GetGridSpacing (const int grid, double * spacing)
{
  if (grid == 0) {
    spacing[0] = this->_model->spacing[0];
  }

}


void BmiMapper::
GetGridOrigin (const int grid, double *origin)
{
  if (grid == 0) {
    origin[0] = this->_model->origin[0];
  }
}


int BmiMapper::
GetGridRank(const int grid)
{
  if (grid == 0)
    return 1;
  else
    return -1;
}


int BmiMapper::
GetGridSize(const int grid)
{
  if (grid == 0)
    return this->_model->shape[0];
  else if (grid == 1)
    return 1;
  else
    return -1;
}


std::string BmiMapper::
GetGridType(const int grid)
{
  if (grid == 0)
    return "uniform_rectilinear";
  else
    return "";
}


void BmiMapper::
GetGridX(const int grid, double *x)
{
  throw mapper::NotImplemented();
}


void BmiMapper::
GetGridY(const int grid, double *y)
{
  throw mapper::NotImplemented();
}


void BmiMapper::
GetGridZ(const int grid, double *z)
{
  throw mapper::NotImplemented();
}


int BmiMapper::
GetGridNodeCount(const int grid)
{
  throw mapper::NotImplemented();
}


int BmiMapper::
GetGridEdgeCount(const int grid)
{
  throw mapper::NotImplemented();
}


int BmiMapper::
GetGridFaceCount(const int grid)
{
  throw mapper::NotImplemented();
}


void BmiMapper::
GetGridEdgeNodes(const int grid, int *edge_nodes)
{
  throw mapper::NotImplemented();
}


void BmiMapper::
GetGridFaceEdges(const int grid, int *face_edges)
{
  throw mapper::NotImplemented();
}


void BmiMapper::
GetGridFaceNodes(const int grid, int *face_nodes)
{
  throw mapper::NotImplemented();
}


void BmiMapper::
GetGridNodesPerFace(const int grid, int *nodes_per_face)
{
  throw mapper::NotImplemented();
}


void BmiMapper::
GetValue (std::string name, void *dest)
{
  void * src = NULL;
  int nbytes = 0;

  src = this->GetValuePtr(name);
  nbytes = this->GetVarNbytes(name);
  
  memcpy (dest, src, nbytes);
}


void *BmiMapper::
GetValuePtr (std::string name)
{
  if (name.compare("TWI") == 0)
    return (void*)this->_model->TWI;
  else if (name.compare("global_deficit") == 0)
    return (void*)(&this->_model->cat_global_deficit);
  else  if (name.compare("porosity") == 0)
    return (void*)(&this->_model->maxsmc);
  else if (name.compare("depth") == 0)
    return (void*)(&this->_model->depth);
  else if (name.compare("grid_gid") == 0)
    return (void*)this->_model->grid_id;
  else if (name.compare("grid_gid_unique") == 0)
    return (void*)this->_model->grid_id_unique;
  else if (name.compare("area_fraction") == 0)
    return (void*)this->_model->grid_area_fraction;
  else if (name.compare("grid_soil_moisture") == 0)
    return (void*)this->_model->grid_soil_moisture;
    else {
    std::stringstream errMsg;
    errMsg << "variable "<< name << " does not exist";
    throw std::runtime_error(errMsg.str());
    return NULL;
  }
}


void BmiMapper::
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


void BmiMapper::
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


void BmiMapper::
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


std::string BmiMapper::
GetComponentName()
{
  return "Mapper BMI";
}


int BmiMapper::
GetInputItemCount()
{
  return this->input_var_name_count;
}


int BmiMapper::
GetOutputItemCount()
{
  return this->output_var_name_count;
}


std::vector<std::string> BmiMapper::
GetInputVarNames()
{
  std::vector<std::string> names;

  for (int i=0; i<this->input_var_name_count; i++)
    names.push_back(this->input_var_names[i]);

  return names;
}


std::vector<std::string> BmiMapper::
GetOutputVarNames()
{
  std::vector<std::string> names;

  for (int i=0; i<this->output_var_name_count; i++)
    names.push_back(this->output_var_names[i]);

  return names;
}


double BmiMapper::
GetStartTime () {
  throw mapper::NotImplemented();
}


double BmiMapper::
GetEndTime () {
  throw mapper::NotImplemented();
}


double BmiMapper::
GetCurrentTime () {
  throw mapper::NotImplemented();
}


std::string BmiMapper::
GetTimeUnits() {
  throw mapper::NotImplemented();
}


double BmiMapper::
GetTimeStep () {
  throw mapper::NotImplemented();
}

#endif
