## Configuration Files
Example configuration files are provided here. To run the given examples see the instructions ([here](https://github.com/NOAA-OWP/SoilMoistureProfiles/blob/ajk/doc_update/RUN.md)).

A detailed description of the parameters for model configuration (i.e., initialize/setup) is provided below. The asterisk (*) denotes calibratable parameters (i.e., `smcmax, b, and satpsi` can be calibrated)

| Variable _____________________ | Datatype _______________ |  Limits _________________ | Units ____ | Role ___ |  Description _____________________________________________________________________|
| :-------- | :-------- | :------ | :----- | :---- |  :----------------------- |
| *smcmax | double  (1D array) | - | - | - | the maximum moisture content (i.e., porosity). Note porosity for layered-based models vary by layers |
| *b | double | - | - | - | the pore size distribution, beta exponent in Clapp-Hornberger function |
| *satpsi | double | - | - | - | saturated capillary head (saturated moisture potential) |
| soil_z | double (1D array) | - | m | - | vertical resolution of the soil moisture profile (depths from the surface) |
| soil_storage_model_depth | double | - | m | - | depth of the soil reservoir model (e.g., CFE). Note: this depth can be different from the depth of the soil moisture profile which is based on `soil_z` |
| soil_storage_model | string | conceptual or layered or topmodel | - | - | if `conceptual`, conceptual models are used for computing the soil moisture profile (e.g., CFE). If `layered`, layered-based soil moisture models are used (e.g., LGAR). If `topmodel`, topmodel's variables are used
| soil_moisture_profile_option | string | constant or linear | - | - | `constant` for layered-constant profile. `linear`  for linearly interpolated values between two consecutive layers. Needed if `soil_storage_model = layered`.
| soil_depth_layers | double (1D array) | - | - | - | Absolute depth of soil layers. Needed if `soil_storage_model = layered`.
| soil_moisture_fraction_depth | double | (0, domain_depth] | m | - | **user specified depth for the soil moisture fraction (default is 40 cm)
| water_table_based_method | string | flux-based or deficit-based | - | - | Needed if `soil_storage_model = topmodel`. `flux-based` uses an iterative scheme, and `deficit-based` uses catchment deficit to compute soil moisture profile

**Note: SoilMoistureProfiles has a bmi output called soil_moisture_fraction (water in the topsoil over total water in the soil column) needs the depth of the topsoil (`soil_moisture_fraction_depth`, default is set to 40 cm -- the top two cells in the soil column of the NWM 3.0)
_________________________________________________________________
