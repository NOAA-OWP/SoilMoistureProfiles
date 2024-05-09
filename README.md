# SoilMoistureProfiles
The soil moisture profiles schemes provide soil moisture distributed over a one-dimensional vertical column and depth to water table. These schemes facilitate coupling among hydrological and thermal models such as (CFE and SFT or LASAM and SFT). It can also be coupled with conceptual hydrologic models to compute rootzone-based actual evapotranspiration. The following three schemes are provided here to compute vertical soil moisture profiles and water table depth.
 * Scheme for conceptual soil reservoirs (e.g., **[CFE](https://github.com/NOAA-OWP/cfe))**
 * Schemes for layered soil reservoirs (e.g., **[LGAR](https://github.com/NOAA-OWP/LGAR-C))**
 * Schemes for topmodel (details are provided below; **[TopModel](https://github.com/NOAA-OWP/topmodel))**

## Conceptual soil reservoir
For conceptual reservoirs, see the following schematic and algorithm. We use the Clap-Hornberger soil moisture characteristic function here, and  **soil_storage** and **soil_storage_chage** are the two inputs passed through BMI.

  ![smp_schematic](https://user-images.githubusercontent.com/15165757/164322224-479477d7-2275-4ce3-a00b-9270cc0d3201.png)

## Layered soil reservoir
For layered soil reservoirs, the two options include
  * constant by layer, and Clap-Horngerger soil moisture characteristic function for the profile below the depth of the last layer
  * linearly interpolated profile between consecutive layers, and Clap-Horngerger soil moisture characteristic function for the profile below the depth of the last layer

## Topmodel based soil reservoir
  * (flux-based method) A method using an iterative scheme to first compute watertable depth and then soil moisture profile ([Blazkova et al. (2002)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2001WR000912))
  * (deficit-based method) A method using catchment deficit to first compute watertable depth and then soil moisture profile ([Franchini et al. (1996)](https://www.sciencedirect.com/science/article/abs/pii/S0022169496800151)


## Build and Run Instructions
Detailed instructions on how to build and run SoilMoistureProfiels (SMP) can be found here [INSTALL](https://github.com/NOAA-OWP/SoilMoistureProfiles/blob/main/INSTALL.md).
  - Test examples highlights
    - Unittest: (see [tests](https://github.com/NOAA-OWP/SoilMoistureProfiles/blob/main/tests/README.md))
    - With CFE: An example coupling CFE with SMP (Soil Moisture Profiles) to compute `watertable` and `soil_moisture_profile`  (see [build/run](https://github.com/NOAA-OWP/SoilMoistureProfiles/blob/main/INSTALL.mdd#pseudo-framework-example-1))
    - With topmodel: An example coupling TopModel to SMP (Soil Moisture Profiles) to compute `watertable` and `soil_moisture_profile` (see [build/run](https://github.com/NOAA-OWP/SoilMoistureProfiles/blob/main/INSTALL.md#pseudo-framework-example-2))
    - Nextgen examples: Realization files for the two above examples (Standalone and with topmodel) are provided [here](https://github.com/NOAA-OWP/SoilMoistureProfiles/blob/main/realizations) (see [build/run](https://github.com/NOAA-OWP/SoilMoistureProfiles/blob/main/INSTALL.md#nextgen-framework-example)).

## Model Configuration File
Detailed description of the parameters for model configuration is provided ([here](https://github.com/NOAA-OWP/SoilMoistureProfiles/tree/main/configs/README.md))

## Getting help
For questions, please contact Ahmad Jan (ahmad.jan(at)noaa.gov), the main developer/maintainer of the repository.

## Known issues or raise an issue
We are constantly looking to improve the model and/or fix bugs as they arise. Please see the Git Issues for known issues or if you want to suggest adding a capability or to report a bug, please open an issue.

## Getting involved
See general instructions to contribute to the model development ([instructions](https://github.com/NOAA-OWP/SoilMoistureProfiles/blob/main/CONTRIBUTING.md)) or simply fork the repository and submit a pull request.
