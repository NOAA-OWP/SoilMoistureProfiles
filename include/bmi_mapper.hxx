#ifndef BMI_MAPPER_H_INCLUDED
#define BMI_MAPPER_H_INCLUDED

/*
  @author: Ahmad Jan (ahmad.jan@noaa.gov)
*/

#include <string.h>
#include "../bmi/bmi.hxx"
#include "soil_moisture_mapping.hxx"

using namespace std;

namespace mapper {
class NotImplemented : public std::logic_error {
  public:
  NotImplemented() : std::logic_error("Not Implemented") { };
};

}

class BmiMapper : public bmixx::Bmi {
  
public:
  BmiMapper() {
    this->input_var_names[0] = "TWI";
    this->input_var_names[1] = "global_deficit";
    this->input_var_names[2] = "porosity";
    this->input_var_names[3] = "depth";
    this->input_var_names[4] = "grid_gid";
    this->input_var_names[5] = "area_fraction";
    
    this->output_var_names[0] = "grid_SMC";
    this->output_var_names[0] = "grid_id_unique";
  };
  
  void Initialize(std::string config_file);
  
  void Update();
  void UpdateUntil(double time);
  void Finalize();
  
  std::string GetComponentName();
  int GetInputItemCount();
  int GetOutputItemCount();
  std::vector<std::string> GetInputVarNames();
  std::vector<std::string> GetOutputVarNames();

  int GetVarGrid(std::string name);
  std::string GetVarType(std::string name);
  int GetVarItemsize(std::string name);
  std::string GetVarUnits(std::string name);
  int GetVarNbytes(std::string name);
  std::string GetVarLocation(std::string name);
  
  double GetCurrentTime();
  double GetStartTime();
  double GetEndTime();
  std::string GetTimeUnits();
  double GetTimeStep();
  
    void GetValue(std::string name, void *dest);
  void *GetValuePtr(std::string name);
  void GetValueAtIndices(std::string name, void *dest, int *inds, int count);
  
  void SetValue(std::string name, void *src);
  void SetValueAtIndices(std::string name, int *inds, int len, void *src);
  
  int GetGridRank(const int grid);
  int GetGridSize(const int grid);
  std::string GetGridType(const int grid);
  
  void GetGridShape(const int grid, int *shape);
  void GetGridSpacing(const int grid, double *spacing);
  void GetGridOrigin(const int grid, double *origin);
  
  void GetGridX(const int grid, double *x);
  void GetGridY(const int grid, double *y);
  void GetGridZ(const int grid, double *z);
  
  int GetGridNodeCount(const int grid);
  int GetGridEdgeCount(const int grid);
  int GetGridFaceCount(const int grid);
  
  void GetGridEdgeNodes(const int grid, int *edge_nodes);
  void GetGridFaceEdges(const int grid, int *face_edges);
  void GetGridFaceNodes(const int grid, int *face_nodes);
  void GetGridNodesPerFace(const int grid, int *nodes_per_face);
private:
  smc_mapping::SoilMoistureMapping _model;
  static const int input_var_name_count = 6;
  static const int output_var_name_count = 2;
  
  std::string input_var_names[6];
  std::string output_var_names[2];
};

#endif
