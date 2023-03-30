#ifndef BMI_SMP_H_INCLUDED
#define BMI_SMP_H_INCLUDED

using namespace std;

#include <string.h>
#include "../bmi/bmi.hxx"
#include "soil_moisture_profile.hxx"


namespace coupler {
  class NotImplemented : public std::logic_error {
  public:
    NotImplemented() : std::logic_error("Not Implemented Function in SoilMoistureProfile") { };
  };

}

class BmiSoilMoistureProfile : public bmixx::Bmi {
public:
  BmiSoilMoistureProfile() {
    this->input_var_names[0] = "soil_storage";
    this->input_var_names[1] = "soil_storage_change";
    this->input_var_names[2] = "soil_moisture_layered";
    this->input_var_names[3] = "soil_depths_layered";
    this->input_var_names[4] = "num_cells_layered";
    this->input_var_names[5] = "Qb_topmodel";        // baseflow in the topmodel
    this->input_var_names[6] = "Qv_topmodel";        // recharge rate to the saturated zone to the un saturated zone
                                                     // in the topmodel
    this->input_var_names[7] = "global_deficit";     // global soil deficit in the topmodel
    
    this->output_var_names[0] = "soil_moisture_profile";  // entire profile of the soil column (1D array)
    this->output_var_names[1] = "soil_water_table";       // depth of the water table from the surface in meters
    this->output_var_names[2] = "soil_moisture_fraction"; // fraction of soil moisture, top 0.4 m (or user specified depth)
  };
  
  void Initialize(std::string config_file);
  
  void Update();
  void UpdateUntil(double time);
  void Finalize();
  void PrintSoilMoistureProfile();
  
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
  soil_moisture_profile::soil_profile_parameters* state;
  static const int input_var_name_count  = 8;
  static const int output_var_name_count = 3;
  
  std::string input_var_names[input_var_name_count];
  std::string output_var_names[output_var_name_count];
  std::string verbosity;
};


#ifdef NGEN
extern "C"
{

  /**
   * Construct this BMI instance as a normal C++ object, to be returned to the framework.
   *
   * @return A pointer to the newly allocated instance.
   */
  BmiSoilMoistureProfile *bmi_model_create() {
    return new BmiSoilMoistureProfile();
  }
  
  /**
   * @brief Destroy/free an instance created with @see bmi_model_create
   * 
   * @param ptr 
   */
  void bmi_model_destroy(BmiSoilMoistureProfile *ptr) {
    delete ptr;
  }

}

#endif

#endif
