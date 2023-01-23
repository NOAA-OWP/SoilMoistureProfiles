#ifndef BMI_MODELCXX_C_INCLUDED
#define BMI_MODELCXX_C_INCLUDED


#include <stdio.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "../bmi/bmi.hxx"
#include "../include/bmi_model_cxx.hxx"
#include "../include/model_cxx.hxx"


/*
  check out note #2, #3, and #4 to see main differences between BMIs for C and C++ based models.
*/

string verbosity="none";

void BmiModelCXX::
Initialize (std::string config_file)
{
  if (config_file.compare("") != 0 )
    this->state = new model_cxx::ModelCXX(config_file); /* note #2/#3: create ModelCXX object with "new".
							   "this->state" is declared in the bmi_model_cxx.hxx
							   model_cxx is a namespace,
							   "::" is scope resolution operator
							   compare it with note #2 and #3 in the bmi_model_c.cxx
							*/
  
}


void BmiModelCXX::
Update()
{
  this->state->Advance(); // note #4: call to Advance method implemented in the model_cxx.cxx to advance the model (with timestep dt)
                         // compare it with note #4 in the bmi_model_c.cxx
}


void BmiModelCXX::
UpdateUntil(double t)
{
  throw coupler::NotImplemented();
}


void BmiModelCXX::
Finalize()
{
  if (this->state)
    this->state->~ModelCXX();
}


int BmiModelCXX::
GetVarGrid(std::string name)
{
  if (name.compare("coefficient_a") == 0 || name.compare("coefficient_b") == 0
      || name.compare("coefficient_c") == 0 || name.compare("coefficient_d") == 0)
    return 0;
  else if (name.compare("root") == 0)
    return 0;
  else
    return -1;
}


std::string BmiModelCXX::
GetVarType(std::string name)
{
  int var_grid = GetVarGrid(name);

  if (var_grid == 0)
    return "double";
  else
    return "";
}


int BmiModelCXX::
GetVarItemsize(std::string name)
{
  std::string var_type = GetVarType(name);

  if (var_type.compare("int") == 0)
    return sizeof(int);
  else if (var_type.compare("double") == 0)
    return sizeof(double);
  else
    return 0;
}


std::string BmiModelCXX::
GetVarUnits(std::string name)
{
  if (name.compare("coefficient_a") == 0 || name.compare("coefficient_b") == 0
      || name.compare("coefficient_c") == 0 || name.compare("coefficient_d") == 0)
    return "none";
  else if (name.compare("root") == 0)
    return "none";
  else
    return "none";
}


int BmiModelCXX::
GetVarNbytes(std::string name)
{
  int itemsize;
  int gridsize;

  itemsize = this->GetVarItemsize(name);
  gridsize = this->GetGridSize(this->GetVarGrid(name));
  return itemsize * gridsize;
}


std::string BmiModelCXX::
GetVarLocation(std::string name)
{
  int var_grid = GetVarGrid(name);

  if (var_grid == 0)
    return "node";
  else
    return "none";
}


void BmiModelCXX::
GetGridShape(const int grid, int *shape)
{
  /*
  if (grid == 2) {
    shape[0] = this->state->shape[0];
  }
  else if (grid == 3) {
    shape[0] =this->state->shape[1];
  }
  */
}


void BmiModelCXX::
GetGridSpacing (const int grid, double * spacing)
{
  /*
  if (grid == 0) {
    spacing[0] =this->state->spacing[0];
  }
  */
}


void BmiModelCXX::
GetGridOrigin (const int grid, double *origin)
{
  /*
  if (grid == 0) {
    origin[0] = 0; this->state->origin[0];
  }
  */
}


int BmiModelCXX::
GetGridRank(const int grid)
{
  if (grid <= 3)
    return 1;
  else
    return -1;
}


int BmiModelCXX::
GetGridSize(const int grid)
{
  if (grid == 0 || grid == 1)
    return 1;
  else
    return -1;
}


std::string BmiModelCXX::
GetGridType(const int grid)
{
  if (grid == 0)
    return "uniform_rectilinear";
  else
    return "";
}


void BmiModelCXX::
GetGridX(const int grid, double *x)
{
  throw coupler::NotImplemented();
}


void BmiModelCXX::
GetGridY(const int grid, double *y)
{
  throw coupler::NotImplemented();
}


void BmiModelCXX::
GetGridZ(const int grid, double *z)
{
  throw coupler::NotImplemented();
}


int BmiModelCXX::
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


int BmiModelCXX::
GetGridEdgeCount(const int grid)
{
  throw coupler::NotImplemented();
}


int BmiModelCXX::
GetGridFaceCount(const int grid)
{
  throw coupler::NotImplemented();
}


void BmiModelCXX::
GetGridEdgeNodes(const int grid, int *edge_nodes)
{
  throw coupler::NotImplemented();
}


void BmiModelCXX::
GetGridFaceEdges(const int grid, int *face_edges)
{
  throw coupler::NotImplemented();
}


void BmiModelCXX::
GetGridFaceNodes(const int grid, int *face_nodes)
{
  throw coupler::NotImplemented();
}


void BmiModelCXX::
GetGridNodesPerFace(const int grid, int *nodes_per_face)
{
  throw coupler::NotImplemented();
}


void BmiModelCXX::
GetValue (std::string name, void *dest)
{
  void * src = NULL;
  int nbytes = 0;

  src = this->GetValuePtr(name);
  nbytes = this->GetVarNbytes(name);
  memcpy (dest, src, nbytes);
}


void *BmiModelCXX::
GetValuePtr (std::string name)
{
  if (name.compare("coefficient_a") == 0)
    return (void*)(&this->state->a);
  else if (name.compare("coefficient_b") == 0)
    return (void*)(&this->state->b);
  else  if (name.compare("coefficient_c") == 0)
    return (void*)(&this->state->c);
  else if (name.compare("coefficient_d") == 0)
    return (void*)(&this->state->d);
  else if (name.compare("root") == 0)
    return (void*)(&this->state->root);
  else {
    std::stringstream errMsg;
    errMsg << "variable "<< name << " does not exist";
    throw std::runtime_error(errMsg.str());
    return NULL;
  }
}


void BmiModelCXX::
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


void BmiModelCXX::
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


void BmiModelCXX::
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


std::string BmiModelCXX::
GetComponentName()
{
  return "CXX-based Model BMI";
}


int BmiModelCXX::
GetInputItemCount()
{
  return this->input_var_name_count;
}


int BmiModelCXX::
GetOutputItemCount()
{
  return this->output_var_name_count;
}


std::vector<std::string> BmiModelCXX::
GetInputVarNames()
{
  std::vector<std::string> names;
  
  for (int i=0; i<this->input_var_name_count; i++)
    names.push_back(this->input_var_names[i]);
  
  return names;
}


std::vector<std::string> BmiModelCXX::
GetOutputVarNames()
{
  std::vector<std::string> names;

  for (int i=0; i<this->output_var_name_count; i++)
    names.push_back(this->output_var_names[i]);

  return names;
}


double BmiModelCXX::
GetStartTime () {
  return 0.0;
}


double BmiModelCXX::
GetEndTime () {
  return 0.0;
}


double BmiModelCXX::
GetCurrentTime () {
  return 0.0;
}


std::string BmiModelCXX::
GetTimeUnits() {
  return "s";
}


double BmiModelCXX::
GetTimeStep () {
  return 0;
}

#endif
