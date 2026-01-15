#ifndef BMI_SMP_H_INCLUDED
#define BMI_SMP_H_INCLUDED

using namespace std;

#include <string.h>
#include "../bmi/bmi.hxx"
#include "soil_moisture_profile.hxx"
#include "vecbuf.hxx"

#include <boost/serialization/access.hpp>

namespace coupler {
  class NotImplemented : public std::logic_error {
  public:
    NotImplemented() : std::logic_error("Not Implemented Function in SoilMoistureProfile") { };
  };

}

class BmiSoilMoistureProfile : public bmi::Bmi {
public:
  BmiSoilMoistureProfile() : m_serialized{}, m_serialized_length(0) {
    this->input_var_names[0]  = "soil_storage";
    this->input_var_names[1]  = "soil_storage_change";
    this->input_var_names[2]  = "num_wetting_fronts";
    this->input_var_names[3]  = "soil_moisture_wetting_fronts";
    this->input_var_names[4]  = "soil_depth_wetting_fronts";
    this->input_var_names[5]  = "Qb_topmodel";        // baseflow in the topmodel
    this->input_var_names[6]  = "Qv_topmodel";        // recharge rate to the saturated zone to the un saturated zone
                                                     // in the topmodel
    this->input_var_names[7]  = "global_deficit";    // global soil deficit in the topmodel

    this->output_var_names[0] = "soil_moisture_profile";  // entire profile of the soil column (1D array)
    this->output_var_names[1] = "soil_water_table";       // depth of the water table from the surface in meters
    this->output_var_names[2] = "soil_moisture_fraction"; // fraction of soil moisture, top 0.4 m (or user specified depth)

    // add calibratable parameters
    this->calib_var_names[0]  = "smcmax";
    this->calib_var_names[1]  = "b";
    this->calib_var_names[2]  = "satpsi";
  };

  void Initialize(std::string config_file) override;

  void Update() override;
  void UpdateUntil(double time) override;
  void Finalize() override;
  void PrintSoilMoistureProfile();

  std::string GetComponentName() override;
  int GetInputItemCount() override;
  int GetOutputItemCount() override;
  std::vector<std::string> GetInputVarNames() override;
  std::vector<std::string> GetOutputVarNames() override;

  int GetVarGrid(std::string name) override;
  std::string GetVarType(std::string name) override;
  int GetVarItemsize(std::string name) override;
  std::string GetVarUnits(std::string name) override;
  int GetVarNbytes(std::string name) override;
  std::string GetVarLocation(std::string name) override;

  double GetCurrentTime() override;
  double GetStartTime() override;
  double GetEndTime() override;
  std::string GetTimeUnits() override;
  double GetTimeStep() override;

  void GetValue(std::string name, void *dest) override;
  void *GetValuePtr(std::string name) override;
  void GetValueAtIndices(std::string name, void *dest, int *inds, int count) override;

  void SetValue(std::string name, void *src) override;
  void SetValueAtIndices(std::string name, int *inds, int len, void *src) override;

  int GetGridRank(const int grid) override;
  int GetGridSize(const int grid) override;
  std::string GetGridType(const int grid) override;

  void GetGridShape(const int grid, int *shape) override;
  void GetGridSpacing(const int grid, double *spacing) override;
  void GetGridOrigin(const int grid, double *origin) override;

  void GetGridX(const int grid, double *x) override;
  void GetGridY(const int grid, double *y) override;
  void GetGridZ(const int grid, double *z) override;

  int GetGridNodeCount(const int grid) override;
  int GetGridEdgeCount(const int grid) override;
  int GetGridFaceCount(const int grid) override;

  void GetGridEdgeNodes(const int grid, int *edge_nodes) override;
  void GetGridFaceEdges(const int grid, int *face_edges) override;
  void GetGridFaceNodes(const int grid, int *face_nodes) override;
  void GetGridNodesPerFace(const int grid, int *nodes_per_face) override;
  void ResetSize (std::string name);

private:
  friend class boost::serialization::access;
  soil_moisture_profile::soil_profile_parameters* state{};
  static const int input_var_name_count  = 8;
  static const int output_var_name_count = 3;
  static const int calib_var_name_count  = 3;

  std::string input_var_names[input_var_name_count];
  std::string output_var_names[output_var_name_count];
  std::string calib_var_names[calib_var_name_count];
  std::string verbosity;

  template<class Archive>
  void serialize(Archive& ar, unsigned int version);
  vecbuf<char> m_serialized;
  uint64_t m_serialized_length;
  void new_serialized();
  void load_serialized(const char* data);
  void free_serialized();
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
