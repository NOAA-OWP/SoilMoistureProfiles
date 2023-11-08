############################################################################################
# Author  : Ahmad Jan
# Contact : ahmad.jan@noaa.gov
# Date    : October 12, 2023
############################################################################################

"""
The script modifies the given realization file and writes a new file with JinjaBMI and SLoTH
module blocks. These modules are mainly used for unit conversion of the variables for parity
with NWM 3.0
"""

import os, sys
import json
import argparse

# mapping CFE output variables to LASAM output variables
cfe_lasam_mapping = {
    "SOIL_STORAGE" : "soil_storage",
    "DIRECT_RUNOFF" : "surface_runoff",
    "POTENTIAL_ET": "potential_evapotranspiration",
    "ACTUAL_ET_mm": "actual_evapotranspiration_mm",
    "Q_OUT": "total_discharge",
    "RAIN_RATE": "precipitation"
    }
    
def main(infile, outfile, input_dir, model_option):
    
    with open(infile) as infile:
        df = json.load(infile)
        
    output_vars_old = df['global']['formulations'][0]['params']["output_variables"]

    # list of variables identified for baseline simulations
    output_vars_lst = [
        "nwm_ponded_depth",
        "soil_water_table",
        "deep_gw_to_channel_flux_output",
        "soil_to_gw_flux_output",
        "giuh_runoff_output",
        "accumulated_evapotranspiration_output",
        "sloth_SNOWT_AVG",
        "actual_evaporation_output",
        "soil_ice_fraction",
        "soil_moisture_fraction",
        "sloth_ACSNOM",
        "sloth_ISNOW",
        "sloth_QRAIN",
        "sloth_FSNO",
        "sloth_SNOWH",
        "sloth_SNLIQ",
        "SNEQV",
        "sloth_QSNOW",
        "TGS",
        "volumetric_soil_moisture_output",
        "sloth_ACCECAN",
        "accumulated_evaporation_output",
        "accumulated_transpiration_output",
        "ground_heat_flux",
        "SOIL_STORAGE",
        "DIRECT_RUNOFF",
        "ACTUAL_ET_mm",
        "POTENTIAL_ET",
        "Q_OUT",
        "RAIN_RATE",
        "land_surface_wind__x_component_of_velocity",
        "land_surface_wind__y_component_of_velocity",
        "land_surface_air__temperature",
        "atmosphere_air_water~vapor__relative_saturation",
        "land_surface_radiation~incoming~longwave__energy_flux",
        "land_surface_radiation~incoming~shortwave__energy_flux",
        "atmosphere_water__liquid_equivalent_precipitation_rate",
        "land_surface_air__pressure"
    ]

    # updated the above list if baseline contains LASAM model (note: officially, LASAM is not part of OWP baseline simulation)
    # essentially it the same simulation only LASAM replacing CFE
    if ("lasam" in model_option):
        for var in output_vars_lst:
            if var in ["SOIL_STORAGE", "DIRECT_RUNOFF", "ACTUAL_ET_mm", "POTENTIAL_ET", "Q_OUT", "RAIN_RATE"]:
                idx = output_vars_lst.index(var)
                output_vars_lst[idx] = cfe_lasam_mapping[var]
            
    df['global']['formulations'][0]['params']["output_variables"] = output_vars_lst
    
    output_vars_new = df['global']['formulations'][0]['params']["output_variables"]
    
    
    output_header_lst = [
        "nwm_ponded_depth[mm]",
        "soil_water_table[m]",
        "deep_gw_to_channel_flux[m3/s]",
        "bottom_flux[m3]",
        "giuh_runoff[m3/s]",
        "AccumEvapotranspiration[mm]",
        "sloth_avg_snow_age",
        "Evaporation[kg m-2 s-1]",
        "soil_ice_fraction[-]",
        "soil_moisture_fraction[-]",
        "sloth_ACSNOM",
        "sloth_ISNOW",
        "sloth_QRAIN",
        "sloth_FSNO",
        "sloth_SNOWH",
        "sloth_SNLIQ",
        "snow_water_equivalent[kg m-2]",
        "sloth_QSNOW",
        "soil_ground_temperature[K]",
        "volumetric_soil_moisture[-]",
        "sloth_ACCECAN",
        "AccumEvaporation[mm]",
        "AccumTranspiration[mm]",
        "ground_heat_flux[W m-2]",
        "soil_storage[m]",
        "direct_runoff[m/h]",
        "AET[m/h]",
        "PET[m/h]",
        "q_out[m/h]",
        "rain_rate[m/h]",
        "U2D[m/s]",
        "V2D[m/s]",
        "T2D[K]",
        "Q2D[kg/kg]",
        "LWDOWN[W m-2]",
        "SWDOWN[W m-2]",
        "RAINRATE[mm s-1]",
        "PSFC[Pa]"
    ]

    # update the global output_header_fields that are written to the output file
    df['global']['formulations'][0]['params']["output_header_fields"] = output_header_lst

    ########################################################################################
    # modify modules in the realization to add SLoTH and JinjaBMI modules
    modules = df['global']['formulations'][0]['params']["modules"]

    # 0 points to sloth model (the vary first model in the row)
    sloth = modules[0]['params']['model_params']

    sloth_exe = modules[0]['params']['library_file']
    
    sloth_model_params = {
        "extinction_coefficient(1,double,1,node)": 0.5,
        "leaf_area_index_input(1,double,1,node)": 10.0,
        "sloth_SNOWT_AVG(1,double,1,node)" : 999.0,
        "sloth_ACSNOM(1,double,1,node)" : 999.0,
        "sloth_ISNOW(1,double,1,node)" : 999.0,
        "sloth_QRAIN(1,double,1,node)" : 999.0,
        "sloth_FSNO(1,double,1,node)" : 999.0,
        "sloth_SNOWH(1,double,1,node)" : 999.0,
        "sloth_SNLIQ(1,double,1,node)" : 999.0,
        "sloth_QSNOW(1,double,1,node)" : 999.0,
        "sloth_ACCECAN(1,double,1,node)" : 999.0,
        "sloth_catchment_area(1,double,m^2,node)" : 999.0
    }
    
    sloth_model_params_updated = {**sloth,** sloth_model_params} #merge dict or use dict1 | dict2
    sloth = modules[0]['params']['model_params'] = sloth_model_params_updated

    ########################################################################################
    jinja_block = {
        "name": "bmi_python",
        "params": {
            "model_type_name": "jinjabmi",
            "python_type": "jinjabmi.Jinja",
            "init_config": "%s/jinjabmi/baseline_support.yml"%input_dir,
            "allow_exceed_end_time": True,
            "main_output_variable": "actual_ET_input",
            "uses_forcing_file": False,
	    "variables_names_map": "to_be_filled"
        }
    }

    if ("cfe" in model_option):
        vars_names_map = {
	    "actual_ET_input": "ACTUAL_ET",
	    "direct_runoff_input": "DIRECT_RUNOFF",
	    "giuh_runoff_input": "GIUH_RUNOFF",
	    "soil_storage_input": "SOIL_STORAGE",
	    "catchment_area_input": "sloth_catchment_area",
	    "deep_gw_to_channel_flux_input": "DEEP_GW_TO_CHANNEL_FLUX",
	    "soil_to_gw_flux_input": "SOIL_TO_GW_FLUX"
	}
    elif ("lasam" in model_option):
        vars_names_map = {
	    "actual_ET_input": "actual_evapotranspiration",
	    "direct_runoff_input": "surface_runoff",
	    "giuh_runoff_input": "giuh_runoff",
	    "soil_storage_input": "soil_storage",
	    "catchment_area_input": "sloth_catchment_area",
	    "deep_gw_to_channel_flux_input": "groundwater_to_stream_recharge",
	    "soil_to_gw_flux_input": "percolation"
	}
        
    jinja_block['params']["variables_names_map"] = vars_names_map    
    modules.append(jinja_block)

    ########################################################################################
    sloth_unit_conversion_block = {
        "name": "bmi_c++",
        "params": {
            "model_type_name": "bmi_c++_sloth",
            "library_file": "%s"%sloth_exe,
            "init_config": "/dev/null",
            "allow_exceed_end_time": True,
            "main_output_variable": "nwm_ponded_depth",
            "uses_forcing_file": False,
            "model_params": "to_be_filled"
        }
    }

    if ("cfe" in model_option):
        vars_names_map = {
	    "nwm_ponded_depth(1,double,mm,node,nwm_ponded_depth_output)": 0.0,
	    "ACTUAL_ET_mm(1,double,mm,node,ACTUAL_ET)": 0.0
        }
    elif ("lasam" in model_option):
        vars_names_map = {
	    "nwm_ponded_depth(1,double,mm,node,nwm_ponded_depth_output)": 0.0,
	    "actual_evapotranspiration_mm(1,double,mm,node,actual_evapotranspiration)": 0.0
        }

    sloth_unit_conversion_block['params']["model_params"] = vars_names_map
    modules.append(sloth_unit_conversion_block)

    # writing to outfile.json
    with open(outfile, "w") as ofile:
        json.dump(df, ofile, indent=4, separators=(", ", ": "), sort_keys=False)

    print ("Output written to %s"%outfile)

    
if __name__ == "__main__":
    
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument("-i",    dest="infile", type=str, required=True, help="the input file")
        parser.add_argument("-o",    dest="outfile", type=str, required=True, help="the input file")
        parser.add_argument("-idir", dest="input_dir", type=str, required=True, help="the inputs dir")
        parser.add_argument("-m",    dest="models_option", type=str, required=True, help="option for models coupling")
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)
        
    main(args.infile, args.outfile, args.input_dir, args.models_option)
