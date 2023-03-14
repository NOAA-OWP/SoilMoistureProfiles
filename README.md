## SoilMoistureProfiles
 Three schemes are provided here to compute vertical soil moisture profiles.
 * Scheme for conceptual soil reserviors **(e.g., CFE)** 
 * Schemes for layered soil reservoirs **(e.g., LGAR)**
 * Schemes for topmodel (details are provided below)
 
 ### Conceptual soil reservior
 For conceptual reservoirs, see the following schematic and algorithm. We use the Clap-Hornberger soil moisture characteristic function here, and  **soil moisture storage** is the main input passed through a BMI.
   
  ![smp_schematic](https://user-images.githubusercontent.com/15165757/164322224-479477d7-2275-4ce3-a00b-9270cc0d3201.png)
  
 ### Layered soil reservior
 For layered soil reserviors, the two options include 
  * constant by layer, and Clap-Horngerger soil moisture characteristic function for the profile below the depth of the last layer
  * linearly interpolated profile between consecutive layers, and Clap-Horngerger soil moisture characteristic function for the profile below the depth of the last layer
  
 ### Topmodel based soil reservoir
  * (flux-based method) A method using an iterative scheme to first compute watertable depth and then soil moisture profile ([Blazkova et al. (2002)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2001WR000912))
  * (deficit-based method) A method using catchment deficit to first compute watertable depth and then soil moisture profile ([Franchini et al. (1996)](https://www.sciencedirect.com/science/article/abs/pii/S0022169496800151)


 ### Compiling and Running SoilMoistureProfiles
 There are two examples (CFE and Topmodel) for running SoilMoistureProfiles as described below. They assume you have [GCC](https://gcc.gnu.org) and [CMAKE](https://cmake.org/) on your machine. These are simple examples and provide a one-way coupling (model --> SoilMoistureProfiles), does not change the functionality of the model, and are only for demonstration purposes.
 
 1. `Option STANDALONE` : An example of using `soil_storage` (conceptual reservoir; CFE) to compute `soil_moisture_profile`
 2. `Option WITHTOPMODEL` : An example of using topmodel's outputs to compute `watertable` and `soil_moisture_profile`
 ```
 git clone https://github.com/NOAA-OWP/SoilMoistureProfiles.git
 cd SoilMoistureProfiles && mkdir build && cd build
 cmake ../ [-DSTANDALONE=ON,-DWITHTOPMODEL=ON] (pick one option, e.g. `cmake ../ -DSTANDALONE=ON`)
 make && cd ..
 ./run_smp.sh [STANDALONE, WITHTOPMODEL] (pick one option) 
 ```

 ### Coupled mode (ngen and pseudo frameworks):
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
| soil_storage_model | string | conceptual or layered or topmodel | - | - | if `conceptual`, conceptual models are used for computing the soil moisture profile (e.g., CFE). If `layered`, layered-based soil moisture models are used (e.g., LGAR). If `topmodel`, topmodel's variables are used
| soil_moisture_profile_option | string | constant or linear | - | - | Only needed if `soil_storage_model = layered`. `constant` for layered-constant profile. `linear`  for linearly interpolated values between two consecutive layers
| soil_moisture_fraction_depth | double | (0, domain_depth] | m | - | *user specified depth for the soil moisture fraction (default is 40 cm)
| water_table_based_method | string | flux-based or deficit-based | - | - | Only needed if `soil_storage_model = topmodel`. `flux-based` uses an iterative scheme, and `deficit-based` uses catchment deficit to compute soil moisture profile

*Note: SoilMoistureProfiles has a bmi output called soil_moisture_fraction (water in the topsoil over total water in the soil column) needs the depth of the topsoil (`soil_moisture_fraction_depth`, default is set to 40 cm -- the top two cells in the soil column of the NWM 3.0)
_________________________________________________________________
