## SoilMoistureProfiles
 Two models are provided here to compute vertical soil moisture profiles.
 * Model for conceptual soil reserviors **(e.g., CFE)** 
 * Model for layered soil reservoirs **(e.g., LGAR)**
 
 ### Conceptual soil reservior
 For conceptual reservoirs, see the following schematic and algorithm. We use the Clap-Hornberger soil moisture characteristic function here, and  **soil moisture storage** is the main input passed through a BMI.
   
  ![smp_schematic](https://user-images.githubusercontent.com/15165757/164322224-479477d7-2275-4ce3-a00b-9270cc0d3201.png)
  
  ### Layered soil reservior
 For layered soil reserviors, the two options include 
  * constant by layer, and Clap-Horngerger soil moisture characteristic function for the profile below the depth of the last layer
  * linearly interpolated profile between consecutive layers, and Clap-Horngerger soil moisture characteristic function for the profile below the depth of the last layer
 
 #### Standalone run:
  * run [make_run_standalone.sh](https://github.com/NOAA-OWP/SoilMoistureProfiles/blob/main/make_run_standalone.sh)

 #### Coupled mode (ngen and pseudo frameworks):
  * Coupling SoilMoistureProfiles to any module (for instance, CFE or SFT) **must** follow these [instructions](https://github.com/NOAA-OWP/SoilFreezeThaw) for building and running. **Note separate instructions are provided for building/running in the ngen framework on the  SoilFreezeThaw repo.**
  * Follow this example: [couple SMP with SFT](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/src/main_cfe_aorc_pet_ftm.cxx)

 
 
_________________________________________________________________
## Description of the parameters in the config file

| Variable _____________________ | Datatype ___________ |  Limits _________________ | Units ____ | Role ___ |  Description _____________________________________________________________________|
| :-------- | :-------- | :------ | :----- | :---- |  :----------------------- |
| smcmax | double   | - | - | - | the maximum moisture content (i.e., porosity) |
| b | double | - | - | - | the pore size distribution, beta exponent in Clapp-Hornberger function |
| satpsi | double | - | - | - | saturated capillary head (saturated moisture potential) |
| soil_z | double (1D array) | - | m | - | vertical resolution of the soil moisture profile (depths from the surface) |
| soil_storage_model_depth | double | - | m | - | depth of the soil reservoir model (e.g., CFE). Note: this depth can be different from the depth of the soil moisture profile which is based on `soil_z` |
| soil_storage_model | string | conceptual or layered | - | - | if `conceptual`, conceptual models are used for computing the soil moisture profile (e.g., CFE). If `layered`, layered-based soil moisture models are used (e.g., LGAR)
| soil_moisture_profile_option | string | constant or linear | - | - | Only needed if `soil_storage_model = layered`. `constant` for layered-constant profile. `linear`  for linearly interpolated values between two consecutive layers
| soil_moisture_fraction_depth | double | (0, domain_depth] | m | - | user specified depth for the soil moisture fraction (default is 40 cm)

_________________________________________________________________
