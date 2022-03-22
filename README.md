# SoilMoistureProfiles
 The code computes vertical soil moisture profiles for conceptual **(e.g., CFE)** and layered **(e.g., LGAR)** soil revervoirs.
 Inputs include **total soil moisture storage** and change in the **soil moisture storage per timestep**
 For conceptual reservoir, the algorithm has two steps.
  * Find water table location using Netwon iterative scheme
  * Use Clap-Hornberger soil moisture characteristic function to compute vertical profile
 For layered soil reserviors, the two options include 
  * constant by layer, and Clap-Horngerger soil moisture characteristic function for the profile below the depth of the last layer
  * linearly interpolated profile between consecutive layers, and Clap-Horngerger soil moisture characteristic function for the profile below the depth of the last layer
  
## Standalone run:
 * run [make_bmi_coupler.sh](https://github.com/NOAA-OWP/SoilMoistureProfiles/blob/main/make_bmi_coupler.sh) to get the executable (e.g., run_bmi_coupler)
 * run `./run_bmi_coupler test/config.txt`

## Coupled mode:
 * Coupling SoilMoistureProfiles to any module (for instance, CFE or SFT) **must** follow these [instructions](https://github.com/NOAA-OWP/SoilFreezeThaw) for building and running
 * Follow this example: [couple SMP with SFT](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/src/main_cfe_aorc_pet_ftm.cxx)

## Parameters in the config file
* `smcmax` (porosity; maximum moisture content), `b`, `satpsi` (saturated matric potential) are needed for the soil moisture characteristic function
* `Z` is the soil layer thickness (for instance, if SMP is coupled with SFT then this `Z` is the vertical resolution of the SFT model)
* `smc_profile` = `conceptual` if conceptual model(e.g., CFE) or `layered` if layered based soil moisture model is used (e.g., LGAR)
* if `smc_profile` = `layered`, then set `smc_profile_option` = `constant` for layered-constant profile or `linear`  for linearly interpolated profile between layers
* `water_table` depth of the initial water table

_________________________________________________________________
_________________________________________________________________

# SoilMoistureMapping
SMCMapping BMI computes and maps local (subcatchment) soil moisture from catchment to NWM grid (1x1 km). The idea of computing local soil moisture is based on the concept of [TopModel](https://github.com/NOAA-OWP/topmodel), where given the global (mean) watershed soil moisture deficit and TWI, local soil moisture deficit is computed. Here, we use the watershed maximum storage capacity and TWI to compute local soil moisture content which is then mapped onto NWM 1x1 km grid using the spatial mapping information provided in the hydrofabric.

  
## Standalone run (for unit test or compare with other gis or python mapping):
 * run [make_bmi_mapper.sh](https://github.com/NOAA-OWP/SoilMoistureProfiles/blob/main/make_bmi_mapper.sh) to get the executable (e.g., run_bmi_mapper)
 * run `./run_bmi_mapper test/config_mapping.txt`

## Coupled mode:
 * Coupling SoilMoistureMapping module to any other models (for instance, TopModel) **must** follow these [instructions](https://github.com/NOAA-OWP/SoilFreezeThaw) for building and running
 * Follow this example: [couple SMP with SFT](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/src/main_cfe_aorc_pet_ftm.cxx)

## Parameters in the config file
* `spatial_map_file` contains NWM GIDs Catchments GIDs and the fraction area (fraction area here means fraction of the 1x1 km grid contained in the subcatchment)  
    <img width="548" alt="Screen Shot 2022-03-22 at 12 52 45 PM" src="https://user-images.githubusercontent.com/15165757/159533177-6d55f665-e77b-439e-b957-5de2db054cab.png">
    
   * okay so we see area fraction greater than 1, it is because the unit test was performed on a 2 x 2 km grid

* `TWI_file` contain catchment GIDs, TWI and the distribution of area corresponding to TWI histogram (dist_area_TWI should sum to 1) 
    <img width="689" alt="Screen Shot 2022-03-22 at 12 49 53 PM" src="https://user-images.githubusercontent.com/15165757/159532540-75a07561-ca88-4687-8187-8be146dfd55a.png">

* `global_deficit` mean watershed deficit
* `porosity` soil pore space
* `depth` soil depth (need to compute maximum soil storage (= porosity * depth))
* `szm` a decay factor of transmissivity with increase in storage deficit with dimensions of length (`m` in the TopModel)
