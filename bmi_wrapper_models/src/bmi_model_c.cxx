#ifndef BMI_SMP_C_INCLUDED
#define BMI_SMP_C_INCLUDED


#include <stdio.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "../bmi/bmi.hxx"
#include "../include/bmi_model_c.hxx"
#include "../include/model_c.hxx"


/*
  check out note #2, #3, and #4 to see main differences between BMIs for C and C++ based models.
*/


string verbosity="none";

void BmiModelC::
Initialize (std::string config_file)
{
  if (config_file.compare("") != 0 ) {
    this->state = new model_c::parameters; /* note #2: allocate "struct parameters" with "new".
                                            "this->state" is declared in the bmi_model_c.hxx
                                            model_c is a namespace, although the model in written in C
					    but since we are using .cxx extension, so we can take the
					    advantage of namespacing"
					    "::" is scope resolution operator
					   */
    
    model_c::ModelC(config_file, state);   /* note #3: pass the config file and state ("struct parameters")
					      to the function ModelC declared in model_c.hxx
					   */
  }
}


void BmiModelC::
Update()
{
  Advance(this->state); // note #4: call to Advance method in model_c.cxx to advance the model (with timestep dt)
                       // For C++ based models this will be "this->state.Advance()"
}


void BmiModelC::
UpdateUntil(double t)
{
  // SoilMoistureProfileUpdate(model);
}


void BmiModelC::
Finalize()
{
  if (this->state)
    delete state;
}


int BmiModelC::
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


std::string BmiModelC::
GetVarType(std::string name)
{
  int var_grid = GetVarGrid(name);

  if (var_grid == 0)
    return "double";
  else
    return "";
}


int BmiModelC::
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


std::string BmiModelC::
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


int BmiModelC::
GetVarNbytes(std::string name)
{
  int itemsize;
  int gridsize;

  itemsize = this->GetVarItemsize(name);
  gridsize = this->GetGridSize(this->GetVarGrid(name));
  return itemsize * gridsize;
}


std::string BmiModelC::
GetVarLocation(std::string name)
{
  int var_grid = GetVarGrid(name);

  if (var_grid <= 3)
    return "node";
  else
    return "none";
}


void BmiModelC::
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


void BmiModelC::
GetGridSpacing (const int grid, double * spacing)
{
  /*
  if (grid == 0) {
    spacing[0] =this->state->spacing[0];
  }
  */
}


void BmiModelC::
GetGridOrigin (const int grid, double *origin)
{
  /*
  if (grid == 0) {
    origin[0] = 0; this->state->origin[0];
  }
  */
}


int BmiModelC::
GetGridRank(const int grid)
{
  if (grid <= 3)
    return 1;
  else
    return -1;
}


int BmiModelC::
GetGridSize(const int grid)
{
  if (grid == 0 || grid == 1)
    return 1;
  else
    return -1;
}


std::string BmiModelC::
GetGridType(const int grid)
{
  if (grid == 0)
    return "uniform_rectilinear";
  else
    return "";
}


void BmiModelC::
GetGridX(const int grid, double *x)
{
  throw coupler::NotImplemented();
}


void BmiModelC::
GetGridY(const int grid, double *y)
{
  throw coupler::NotImplemented();
}


void BmiModelC::
GetGridZ(const int grid, double *z)
{
  throw coupler::NotImplemented();
}


int BmiModelC::
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


int BmiModelC::
GetGridEdgeCount(const int grid)
{
  throw coupler::NotImplemented();
}


int BmiModelC::
GetGridFaceCount(const int grid)
{
  throw coupler::NotImplemented();
}


void BmiModelC::
GetGridEdgeNodes(const int grid, int *edge_nodes)
{
  throw coupler::NotImplemented();
}


void BmiModelC::
GetGridFaceEdges(const int grid, int *face_edges)
{
  throw coupler::NotImplemented();
}


void BmiModelC::
GetGridFaceNodes(const int grid, int *face_nodes)
{
  throw coupler::NotImplemented();
}


void BmiModelC::
GetGridNodesPerFace(const int grid, int *nodes_per_face)
{
  throw coupler::NotImplemented();
}


void BmiModelC::
GetValue (std::string name, void *dest)
{
  void * src = NULL;
  int nbytes = 0;

  src = this->GetValuePtr(name);
  nbytes = this->GetVarNbytes(name);
  memcpy (dest, src, nbytes);
}


void *BmiModelC::
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


void BmiModelC::
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


void BmiModelC::
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


void BmiModelC::
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


std::string BmiModelC::
GetComponentName()
{
  return "C-based Model BMI";
}


int BmiModelC::
GetInputItemCount()
{
  return this->input_var_name_count;
}


int BmiModelC::
GetOutputItemCount()
{
  return this->output_var_name_count;
}


std::vector<std::string> BmiModelC::
GetInputVarNames()
{
  std::vector<std::string> names;
  
  for (int i=0; i<this->input_var_name_count; i++)
    names.push_back(this->input_var_names[i]);
  
  return names;
}


std::vector<std::string> BmiModelC::
GetOutputVarNames()
{
  std::vector<std::string> names;

  for (int i=0; i<this->output_var_name_count; i++)
    names.push_back(this->output_var_names[i]);

  return names;
}


double BmiModelC::
GetStartTime () {
  return 0.0;
}


double BmiModelC::
GetEndTime () {
  return 0.0;
}


double BmiModelC::
GetCurrentTime () {
  return 0.0;
}


std::string BmiModelC::
GetTimeUnits() {
  return "s";
}


double BmiModelC::
GetTimeStep () {
  return 0;
}

#endif
