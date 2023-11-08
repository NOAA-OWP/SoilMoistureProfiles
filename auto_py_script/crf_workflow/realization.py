############################################################################################
# Author  : Ahmad Jan
# Contact : ahmad.jan@noaa.gov
# Date    : September 28, 2023
############################################################################################

"""
The script generates realization files for a given model coupling option
 - inputs  : see main function for inputs (taken as arguments)
 - outputs : writes nextgen realization file for the basin
"""


import os, sys
import glob
import json
import subprocess
import argparse
import pandas as pd
import geopandas as gpd
import numpy as np

#############################################################################
# module for potential evapotranspiration block in the nextgen realization file 
# @param input_dir : input directory of the PET config files
# @param model_exe : path to PET executable
#############################################################################
def get_pet_block(model_exe, input_dir):
    
    block = {
	"name": "bmi_c",
	"params": {
            "model_type_name": "bmi_c_pet",
            "library_file": model_exe,
            "forcing_file": "",
            "init_config": os.path.join(input_dir, 'pet/{{id}}_config_pet.txt'), 
            "allow_exceed_end_time": "true",
            "main_output_variable": "water_potential_evaporation_flux",
            "registration_function":"register_bmi_pet",
            "uses_forcing_file": "false"
	}
    }

    return block

#############################################################################
# module for NOAH-OWP-Modular (NOM) block in the nextgen realization file 
# @param input_dir : input directory of the NOM config files
# @param model_exe : path to NOM executable
#############################################################################
def get_noah_owp_modular_block(model_exe, input_dir):
    block = {
        "name": "bmi_fortran", 
        "params": {
            "name": "bmi_fortran", 
            "model_type_name": "NoahOWP", 
            "main_output_variable": "QINSUR",
            "library_file": model_exe,
            "init_config": os.path.join(input_dir, 'nom/{{id}}_config_nom.input'),
            "allow_exceed_end_time": True,
            "fixed_time_step": False,
            "uses_forcing_file": False,
            "variables_names_map": {
                "PRCPNONC": "atmosphere_water__liquid_equivalent_precipitation_rate",
                "Q2": "atmosphere_air_water~vapor__relative_saturation",
                "SFCTMP": "land_surface_air__temperature",
                "UU": "land_surface_wind__x_component_of_velocity",
                "VV": "land_surface_wind__y_component_of_velocity",
                "LWDN": "land_surface_radiation~incoming~longwave__energy_flux",
                "SOLDN": "land_surface_radiation~incoming~shortwave__energy_flux",
                "SFCPRS": "land_surface_air__pressure"
            }
        }
    }
    
    return block

#############################################################################
# module for CFE block in the nextgen realization file 
# @param input_dir : input directory of the CFE config files
# @param model_exe : path to CFE executable
# @param cfe_standalone : if true, additional parameters are mapped from the SLoTH BMI
#############################################################################
def get_cfe_block(model_exe, input_dir, cfe_standalone):
    
    block = {
        "name": "bmi_c",
        "params": {
            "name": "bmi_c", 
            "model_type_name": "CFE", 
            "main_output_variable": "Q_OUT",
            "library_file": model_exe,
            "init_config": os.path.join(input_dir, 'cfe/{{id}}_config_cfe.txt'), 
            "allow_exceed_end_time": True,
            "fixed_time_step": False,
            "uses_forcing_file": False,
            "variables_names_map": {
                "atmosphere_water__liquid_equivalent_precipitation_rate": "QINSUR",
                "water_potential_evaporation_flux": "EVAPOTRANS"
            },
            "registration_function": "register_bmi_cfe"
        }
    }

    if (cfe_standalone):
        sub_map = {
            "ice_fraction_schaake": "ice_fraction_schaake",
            "ice_fraction_xinanjiang": "ice_fraction_xinanjiang",
            "soil_moisture_profile": "soil_moisture_profile"
        }
        var_name_map = block["params"]["variables_names_map"]
        var_name_map.update(sub_map)
        block["params"]["variables_names_map"] = var_name_map
        
    return block

#############################################################################
# module for topmodel block in the nextgen realization file 
# @param input_dir : input directory of the PET config files
# @param model_exe : path to PET executable
#############################################################################
def get_topmodel_block(model_exe, input_dir):
    block = {
        "name": "bmi_c",
        "params": {
            "name": "bmi_c", 
            "model_type_name": "TOPMODEL", 
            "main_output_variable": "Qout",
            "library_file": model_exe,
            "init_config": os.path.join(input_dir, 'topmodel/{{id}}_config_topmodel.run'),
            "allow_exceed_end_time": True,
            "fixed_time_step": False,
            "uses_forcing_file": False,
            "variables_names_map": {
                "atmosphere_water__liquid_equivalent_precipitation_rate": "QINSUR",
                "water_potential_evaporation_flux": "EVAPOTRANS"
            },
            "registration_function": "register_bmi_topmodel"
        }
    }

#############################################################################
# module for Soil Freeze-thaw (SFT) block in the nextgen realization file 
# @param input_dir : input directory of the SFT config files
#############################################################################
def get_sft_block(model_exe, input_dir):
    block = {
        "name": "bmi_c++",
        "params": {
            "name": "bmi_c++",
            "model_type_name": "SFT", 
            "main_output_variable": "num_cells",
            "library_file": model_exe,
            "init_config": os.path.join(input_dir, 'sft/{{id}}_config_sft.txt'),
            "allow_exceed_end_time": True, 
            "uses_forcing_file": False,
            "variables_names_map": {
                "ground_temperature" : "TGS"
            }
        }
    }
    
    return block

#############################################################################
# module for Soil Moisture Profiles (SMP) block in the nextgen realization file 
# @param input_dir : input directory of the SMP config files
# @param model_exe : path to SMP executable
#############################################################################
def get_smp_block(model_exe, input_dir, coupled_models):
    block = {
        "name": "bmi_c++",
        "params": {
            "name": "bmi_c++", 
            "model_type_name": "SMP", 
            "main_output_variable": "soil_water_table",
            "library_file": model_exe,
            "init_config": os.path.join(input_dir, 'smp/{{id}}_config_smp.txt'),
            "allow_exceed_end_time": True,
            "uses_forcing_file": False,
            "variables_names_map": {
                "soil_storage": "SOIL_STORAGE",
		"soil_storage_change": "SOIL_STORAGE_CHANGE"
            }
        }
    }

    name_map = dict()

    if (coupled_models == "nom_cfe_smp_sft" or coupled_models == "cfe_smp" or
        coupled_models == "nom_cfe_smp"):
        name_map = {
            "soil_storage": "SOIL_STORAGE",
	    "soil_storage_change": "SOIL_STORAGE_CHANGE"
        }
    elif (coupled_models == "nom_lasam_smp_sft"):
        name_map = {
            "soil_storage": "sloth_soil_storage",
	    "soil_storage_change": "sloth_soil_storage_change",
	    "soil_moisture_wetting_fronts" : "soil_moisture_wetting_fronts",
	    "soil_depth_wetting_fronts" : "soil_depth_wetting_fronts",
	    "num_wetting_fronts" : "soil_num_wetting_fronts"
        }
    else:
        print ("coupled_models name should be nom_cfe_smp_sft or nom_lasam_smp_sft, provided is ", coupled_models)
        quit()
        
    model_map = block["params"]
    model_map["variables_names_map"] = name_map

    return block

#############################################################################
# module for Lumped Arid/semi-arid model (LASAM) block in the nextgen realization file 
# @param input_dir : input directory of the LASAM config files
# @param model_exe : path to LASAM executable
#############################################################################
def get_lasam_block(model_exe, input_dir):
    block = {
        "name": "bmi_c++",
        "params": {
            "name": "bmi_c++",
            "model_type_name": "LASAM",
            "main_output_variable": "precipitation_rate",
            "library_file": model_exe,
            "init_config": os.path.join(input_dir, 'lasam/{{id}}_config_lasam.txt'),
            "allow_exceed_end_time": True,
            "uses_forcing_file": False,
            "variables_names_map": {
                "precipitation_rate" : "QINSUR",
                "potential_evapotranspiration_rate": "EVAPOTRANS"
            }
        }
    }

    return block

#############################################################################
# module for SLoTH block in the nextgen realization file
# @param model_exe      : path to SLoTH executable
# @param coupled_models : model option (SLoTH provides inputs to several models, depending on the coupling)
#############################################################################
def get_sloth_block(model_exe, coupled_models):
    
    block = {
        "name": "bmi_c++",
        "params": {
            "name": "bmi_c++", 
            "model_type_name": "SLOTH", 
            "main_output_variable": "z", 
            "library_file": model_exe, 
            "init_config": '/dev/null',
            "allow_exceed_end_time": True, 
            "fixed_time_step": False, 
            "uses_forcing_file": False,
        }
    }

    params = dict()
    
    if (coupled_models in ["nom_cfe", "cfe"]):
        params = {
            "ice_fraction_schaake(1,double,m,node)": 0.0,
            "ice_fraction_xinanjiang(1,double,1,node)": 0.0,
	    "soil_moisture_profile(1,double,1,node)": 0.0
        }
    elif (coupled_models == "nom_lasam"):
        params = {
            "soil_temperature_profile(1,double,K,node)" : 275.15
        }
    elif (coupled_models == "nom_cfe_smp_sft"):
        params = {
            "soil_moisture_wetting_fronts(1,double,1,node)": 0.0,
	    "soil_depth_wetting_fronts(1,double,1,node)": 0.0,
	    "num_wetting_fronts(1,int,1,node)": 1.0,
	    "Qb_topmodel(1,double,1,node)": 0.0,
	    "Qv_topmodel(1,double,1,node)": 0.0,
	    "global_deficit(1,double,1,node)": 0.0
        }
    elif (coupled_models == "nom_lasam_smp_sft"):
        params = {
            "sloth_soil_storage(1,double,m,node)" : 1.0E-10,
            "sloth_soil_storage_change(1,double,m,node)" : 0.0,
	    "Qb_topmodel(1,double,1,node)": 0.0,
	    "Qv_topmodel(1,double,1,node)": 0.0,
	    "global_deficit(1,double,1,node)": 0.0
        }
    else:
        msg = "coupled_models name should be nom_cfe, or nom_cfe_smp_sft or nom_lasam_smp_sft, provided is "+ coupled_models
        sys.exit(msg)
    
    model_params = block["params"]
    model_params["model_params"] = params
    
    return block

#############################################################################
# module for JinjaBMI block and unit conversion in the nextgen realization file (for baseline simulations)
# @param model_exe : path to SLoTH executable
# @param input_dir : input directory of the LASAM config files
#############################################################################
def get_jinjabmi_unit_conversion_block(model_exe, input_dir):
    
    block_jinjabmi = {
        "name": "bmi_python",
        "params": {
            "model_type_name": "jinjabmi",
            "python_type": "jinjabmi.Jinja",
            "init_config": os.path.join(input_dir, "jinjabmi/baseline_support.yml"),
            "allow_exceed_end_time": True,
            "main_output_variable": "actual_ET_input",
            "uses_forcing_file": False,
	    "variables_names_map": {
	        "actual_ET_input": "ACTUAL_ET",
	        "direct_runoff_input": "DIRECT_RUNOFF",
	        "giuh_runoff_input": "GIUH_RUNOFF",
	        "soil_storage_input": "SOIL_STORAGE",
	        "catchment_area_input": "sloth_catchment_area",
	        "deep_gw_to_channel_flux_input": "DEEP_GW_TO_CHANNEL_FLUX",
	        "soil_to_gw_flux_input": "SOIL_TO_GW_FLUX",
	        "giuh_runoff_input": "GIUH_RUNOFF"
	    }
        }
    }

    block_unit_conversion = {
	"name": "bmi_c++",
	"params": {
            "model_type_name": "bmi_c++_sloth",
            "library_file": model_exe,
            "init_config": "/dev/null",
            "allow_exceed_end_time": True,
            "main_output_variable": "nwm_ponded_depth",
            "uses_forcing_file": False,
            "model_params": {
		"nwm_ponded_depth(1,double,mm,node,nwm_ponded_depth_output)": 0.0,
		"ACTUAL_ET_mm(1,double,mm,node,ACTUAL_ET)": 0.0
            }
	}
    }

    block = [ block_jinaBMI, block_unit_conversion]

    return block

#############################################################################
# module that calls all module blocks and assembles/writes full realization block/file
# @param ngen_dir        : path to nextgen directory
# @param forcing_dir : forcing data directory containing data for each catchment
# @param input_dir        : input directory (config files of all models exist here under subdirectories)
# @param realization file : name of the output realization file
# @param coupled_models : models coupling option (pre-defined names; see main.py)
# @param runoff_schame  : surface runoff schemes - Options = Schaake or Xinanjiang (For CFE and SFT)
# @param simulation_time  : dictionary containing simulation start/end time
# @param baseline_casae   : boolean (if true, baseline scenario realization file is requested)
#############################################################################
def write_realization_file(ngen_dir, forcing_dir, input_dir, realization_file,
                           coupled_models, runoff_scheme, simulation_time, baseline_case):

    lib_file = {}
    extern_path = os.path.join(ngen_dir, 'extern')
    models = os.listdir(extern_path)
    lib_files = {}
    # , include topmodel later
    for m in models:
        if m in ['SoilFreezeThaw', 'cfe', 'SoilMoistureProfiles', 'LASAM', 'LGAR-C', \
                 'sloth', 'evapotranspiration', 'noah-owp-modular', 'topmodel']:
            path_m = os.path.join(os.path.join(extern_path,m), "cmake_build")
            if (os.path.exists(path_m)):
                exe_m = glob.glob(os.path.join(path_m,'*bmi.dylib'))
                if exe_m == []:
                    exe_m = glob.glob(os.path.join(path_m,'*model.dylib'))
            
                exe_m = exe_m[0].split('extern')
                exe_m = exe_m[1].split('.')
                #lib_files[m] = './ngen/extern'+str(exe_m[0])
                lib_files[m] = os.path.join(ngen_dir, 'extern'+str(exe_m[0]))
            else:
                lib_files[m] = ""

    print ("\n********** Models executables under extern **************")
    for key, value in lib_files.items():
        print ("Model: ", key, ",", value)
    
    # noah
    nom_block = dict()
    if ("nom" in coupled_models):
        assert(lib_files['noah-owp-modular'] != "")
        nom_block = get_noah_owp_modular_block(lib_files['noah-owp-modular'], input_dir)
    

    # cfe
    cfe_block = dict()
    cfe_standalone= False
    if ("cfe" in coupled_models):
        assert (lib_files['cfe'] != "")
        cfe_block = get_cfe_block(lib_files['cfe'], input_dir, cfe_standalone=False)
    
    # topmodel
    topmodel_block = dict()
    if ("topmodel" in coupled_models):
        assert (lib_files['topmodel'] != "")
        topmodel_block = get_topmodel_block(lib_files['topmodel'], input_dir)
    
    # sloth
    sloth_block = dict()
    if ("cfe" in coupled_models or "smp" in coupled_models or "sft" in coupled_models):
        assert (lib_files['sloth'] != "")
        sloth_block = get_sloth_block(lib_files['sloth'], coupled_models)

   
    # sft
    sft_block = dict()
    if ('sft' in coupled_models):
        assert (lib_files['SoilFreezeThaw'] != "")
        sft_block = get_sft_block(lib_files['SoilFreezeThaw'], input_dir=input_dir)
    
    # smp
    smp_block = dict()
    if ('smp' in coupled_models):
        assert (lib_files['SoilMoistureProfiles'] != "")
        smp_block = get_smp_block(lib_files['SoilMoistureProfiles'], input_dir, coupled_models)

    # lasam
    lasam_block = dict()
    if ('lasam' in coupled_models):
        if ("LASAM" in lib_files.keys()):
            lasam_block = get_lasam_block(lib_files['LASAM'], input_dir)
        elif ("LGAR-C" in lib_files.keys()):
            lasam_block = get_lasam_block(lib_files['LGAR-C'], input_dir)

    pet_block = dict()
    if (lib_files['evapotranspiration'] != ""):
        pet_block = get_pet_block(lib_files['evapotranspiration'], input_dir)


    ##########################################################################
    ## now create root block and add time block and modules block to it
    # global configuration
    root = {
        "time": {
            "start_time": simulation_time['start_time'],
            "end_time": simulation_time['end_time'],
            "output_interval": 3600
        },
        "global": {
            "formulations": "to_be_filled_in",
            "forcing": {
                "file_pattern": ".*{{id}}.*.csv",
                "path": forcing_dir,
                "provider": "CsvPerFeature"
            }
        },
        "output_root": "./ngen_output"
    }

    
    # routing block
    routing_block = {
        "routing": {
            "t_route_connection_path": os.path.join(input_dir, "extern/t-route"),
            "t_route_config_file_with_path": os.path.join(input_dir, "troute_config.yaml")            
        }
    }
    #root.update(routing_block)


    
    # fill the formulation block "to_be_filled_in" for the global list
    global_block = {
        "name": "bmi_multi", 
        "params": {
            "name": "bmi_multi",
            "model_type_name": "",
            "init_config": "",
            "allow_exceed_end_time": False,
            "fixed_time_step": False, 
            "uses_forcing_file": False
        }
    }

    model_type_name = ""
    main_output_variable = ""
    modules = []

    if (coupled_models == "cfe"):
        model_type_name = "CFE"
        main_output_variable = "Q_OUT"
        modules = [sloth_block, pet_block, cfe_block]
        output_variables = ["RAIN_RATE", "DIRECT_RUNOFF", "GIUH_RUNOFF", "NASH_LATERAL_RUNOFF",
	                    "DEEP_GW_TO_CHANNEL_FLUX", "Q_OUT", "SOIL_STORAGE", "POTENTIAL_ET", "ACTUAL_ET"]
        output_header_fields = ["rain_rate", "direct_runoff", "giuh_runoff", "nash_lateral_runoff",
                                "deep_gw_to_channel_flux", "q_out", "soil_storage",  "PET", "AET"]
    elif (coupled_models == "nom_cfe"):
        model_type_name = "NOM_CFE"
        main_output_variable = "Q_OUT"
        modules = [sloth_block, nom_block, cfe_block]
        output_variables = ["RAIN_RATE", "DIRECT_RUNOFF", "GIUH_RUNOFF", "NASH_LATERAL_RUNOFF",
	                    "DEEP_GW_TO_CHANNEL_FLUX", "Q_OUT", "SOIL_STORAGE", "POTENTIAL_ET", "ACTUAL_ET"]
        output_header_fields = ["rain_rate", "direct_runoff", "giuh_runoff", "nash_lateral_runoff",
                                "deep_gw_to_channel_flux", "q_out", "soil_storage", "PET", "AET"]
    elif (coupled_models == "nom_lasam"):
        model_type_name = "NOM_LASAM"
        main_output_variable = "total_discharge"
        modules = [sloth_block, nom_block, lasam_block]

        output_variables = ["TGS", "precipitation", "potential_evapotranspiration", "actual_evapotranspiration", 
                            "soil_storage", "surface_runoff", "giuh_runoff", "groundwater_to_stream_recharge",  "percolation", "total_discharge", 
                            "infiltration"] 
        output_header_fields = ["ground_temperature", "rain_rate", "PET_rate", "actual_ET",  
                                "soil_storage", "direct_runoff", "giuh_runoff", "deep_gw_to_channel_flux", "soil_to_gw_flux", "q_out",
                                "infiltration"]
    elif (coupled_models == "nom_topmodel"):
        model_type_name = "NOM_TOPMODEL"
        main_output_variable = "Qout"
        modules = [nom_block, topmodel_block]
        output_variables = ["RAIN_RATE", "Qout", "soil_water__domain_volume_deficit","land_surface_water__runoff_mass_flux"]
        output_header_fields = ["rain_rate", "qout", "soil_deficit", "direct_runoff"]
    elif (coupled_models == "nom_cfe_smp_sft"):
        model_type_name = "NOM_CFE_SMP_SFT"
        main_output_variable = "Q_OUT"
        modules = [sloth_block, nom_block,cfe_block, smp_block, sft_block]
        
        output_variables = ["soil_ice_fraction", "TGS", "RAIN_RATE", "DIRECT_RUNOFF", "GIUH_RUNOFF", "NASH_LATERAL_RUNOFF",
	                    "DEEP_GW_TO_CHANNEL_FLUX", "Q_OUT", "SOIL_STORAGE", "POTENTIAL_ET", "ACTUAL_ET", "soil_moisture_fraction","ice_fraction_schaake"]
        output_header_fields = ["soil_ice_fraction", "ground_temperature", "rain_rate", "direct_runoff", "giuh_runoff", "nash_lateral_runoff",
                                "deep_gw_to_channel_flux", "q_out", "soil_storage", "PET", "AET", "soil_moisture_fraction","ice_fraction_schaake"]
        
        if (runoff_scheme == "Xinanjiang"):
            output_variables[-1] = "ice_fraction_xinanjiang"
            output_header_fields[-1] = "ice_fraction_xinanjiang"
        
    elif (coupled_models == "nom_lasam_smp_sft"):
        model_type_name = "NOM_LASAM_SMP_SFT"
        main_output_variable = "total_discharge"
        modules = [sloth_block, nom_block, lasam_block, smp_block, sft_block]

        output_variables = ["soil_ice_fraction", "TGS", "precipitation", "potential_evapotranspiration", "actual_evapotranspiration", 
                            "soil_storage", "surface_runoff", "giuh_runoff", "groundwater_to_stream_recharge",  "percolation", "total_discharge", 
                            "infiltration", "soil_moisture_fraction"] 
        output_header_fields = ["soil_ice_fraction", "ground_temperature", "rain_rate", "PET_rate", "actual_ET",  
                                "soil_storage", "direct_runoff", "giuh_runoff", "deep_gw_to_channel_flux", "soil_to_gw_flux", "q_out",
                                "infiltration", "soil_moisture_fraction"]

    assert len(output_variables) == len(output_header_fields)
    
    global_block["params"]["model_type_name"] = model_type_name
    global_block["params"]["main_output_variable"] = main_output_variable
    global_block["params"]["output_variables"] = output_variables
    global_block["params"]["output_variables"] = output_variables
    global_block["params"]["output_header_fields"] = output_header_fields
    global_block["params"]["modules"] = modules


    # replace formulations block
    root["global"]["formulations"] = [global_block]
    
    # save realization file as .json
    with open(realization_file, 'w') as outfile:
        print (realization_file)
        json.dump(root, outfile, indent=4, separators=(", ", ": "), sort_keys=False)

    

def main():

    try:
        parser = argparse.ArgumentParser()
        parser.add_argument("-ngen", dest="ngen_dir",   type=str, required=True, help="the giuh files directory")
        parser.add_argument("-f", dest="forcing_dir",   type=str, required=False, help="the forcing files directory")
        parser.add_argument("-i", dest="input_dir",     type=str, required=True, help="the input files directory")
        parser.add_argument("-m", dest="models_option", type=str, required=True, help="option for models coupling")
        parser.add_argument("-r", dest="runoff_scheme", type=str, required=False, help="option for runoff scheme")
        parser.add_argument("-t",    dest="time",       type=json.loads, required=True,  help="simulation start/end time") 
        parser.add_argument("-b", dest="baseline_case", type=str, required=False, help="option for baseline case", default=False) 
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    if (not os.path.exists(args.forcing_dir)):
        sys.exit("Forcing directory does not exist: " + args.forcing_dir)
        
    cfe_dir = os.path.join(args.input_dir,"cfe")
    if ( ('cfe' in args.models_option) and (not os.path.exists(cfe_dir)) ):
        print ("cfe config files directory does not exist under " + args.input_dir)
        sys.exit(0)
    
    sft_file = os.path.join(args.input_dir,"sft")
    if ( ('sft' in args.models_option) and (not os.path.exists(sft_file)) ):
        print ("sft config files directory does not exist under " + args.input_dir)
        sys.exit(0)

    smp_file = os.path.join(args.input_dir,"smp")
    if ( ('smp' in args.models_option) and (not os.path.exists(smp_file)) ):
        print ("smp config files directory does not exist under " + args.input_dir)
        sys.exit(0)

    lasam_file = os.path.join(args.input_dir,"lasam")    
    if ( ("lasam" in args.models_option) and (not os.path.exists(lasam_file)) ):
        print ("lasam config files directory does not exist under " + args.input_dir)
        sys.exit(0)

    write_realization_file(
        ngen_dir     = args.ngen_dir,
        forcing_dir  = args.forcing_dir,
        input_dir    = args.input_dir, 
        realization_file = os.path.join(os.getcwd(), "realization_%s.json"%args.models_option),
        coupled_models   = args.models_option,
        runoff_scheme    = args.runoff_scheme,
        simulation_time  = args.time,
        baseline_case    = args.baseline_case
    )


if __name__ == "__main__":
    main()    
