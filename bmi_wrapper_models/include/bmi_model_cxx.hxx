#ifndef BMI_MODELCXX_H_INCLUDED
#define BMI_MODELCXX_H_INCLUDED

using namespace std;

#include <string.h>
#include "../bmi/bmi.hxx"
#include "model_cxx.hxx"

/*
  check out note #1 to see the main difference between BMIs for C and C++ based models.
*/


namespace coupler {
class NotImplemented : public std::logic_error {
  public:
  NotImplemented() : std::logic_error("Not Implemented Function in Model CXX") { };
};

}

class BmiModelCXX : public bmixx::Bmi {

public:
  BmiModelCXX() {
    this->input_var_names[0] = "coefficient_a";
    this->input_var_names[1] = "coefficient_b";
    this->input_var_names[2] = "coefficient_c";
    this->input_var_names[3] = "coefficient_d";
      
    this->output_var_names[0] = "root";
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
  model_cxx::ModelCXX* state; // note #1: pointer object to the class ModelCXX in model_cxx.hxx
  static const int input_var_name_count = 4;
  static const int output_var_name_count = 1;
  
  std::string input_var_names[input_var_name_count];
  std::string output_var_names[output_var_name_count];
};


#ifdef NGEN
extern "C"
{

  /**
   * Construct this BMI instance as a normal C++ object, to be returned to the framework.
   *
   * @return A pointer to the newly allocated instance.
   */
  BmiModelCXX *bmi_model_create()
  {
    return new BmiModelCXX();
  }
  
  /**
   * @brief Destroy/free an instance created with @see bmi_model_create
   * 
   * @param ptr 
   */
  void bmi_model_destroy(BmiModelCXX *ptr)
  {
    delete ptr;
  }
  
}

#endif

#endif
