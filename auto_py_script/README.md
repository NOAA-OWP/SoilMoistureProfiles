## Instructions for computing Geomorphological Instantaneous Unit Hydrograph (GIUH) and Topographic Wetness Index (TWI)
This set of scripts offers a workflow for establishing basin-scale simulations from the ground up 
and executing them within the NextGen framework.

- R-based scripts, leveraging the [WhiteBox](https://www.whiteboxgeo.com/manual/wbw-user-manual/book/tool_help.html) tool, a geospatial data analysis software,
  are provided for the computation of GIUH and TWI.
  Detailed instruction are available [here](https://github.com/ajkhattak/SoilMoistureProfiles/blob/ajk/auto_py_script/auto_py_script/R/main.R)
- Python-based scripts are provided for generating model(s) (standalone and coupled) configuration files and the
  NextGen realization file. Detailed instruction can be found [here](https://github.com/ajkhattak/SoilMoistureProfiles/tree/ajk/auto_py_script/auto_py_script/crf_workflow/main.py)

### NOTE
This workflow does not download basin's forcing data. Users are required to provide the forcing data. 
General instructions on how to download forcing data are documented [here](https://github.com/ajkhattak/SoilMoistureProfiles/blob/ajk/auto_py_script/auto_py_script/FORCING.md).

### Plotting
To visualize and compare WhiteBox-based vs TauDEM-based GIUH and TWI, follow the steps below. A basin in Maine (gauge ID `1033000`) is provided as an example. This script will generate visualizations comparing the GIUH results obtained from WhiteBox and TauDEM for a few selected catchments within the basin.

#### GIUH Comparison:
Run the script `plot_giuh.R` located in the R directory.

#### TWI Comparison:
Run the script `plot_twi.R` located in the R directory.
