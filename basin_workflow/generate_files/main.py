############################################################################################
# Author  : Ahmad Jan Khattak
# Contact : ahmad.jan@noaa.gov
# Date    : October 11, 2023 
############################################################################################


import os,sys
import subprocess
import pandas as pd
import glob
import shutil
import re
import geopandas as gpd

# Note #1: from the command line just run 'python path_to/main.py'
# Note #2: make sure to adjust the following required arguments
# Note #3: several model coupling options are available, the script currently supports a few of them, for full list see
#          below the coupled_models_option


############################# Required Arguments ###################################
### Specify the following four directories (change according to your settings)
# root_dir     : geopackage(s) directory
# forcing_dir  : lumped forcings directory (pre-downloaed forcing data for each catchment (.csv or .nc))
# ngen_dir     : nextgen directory path
# output_dir   : output directory path (all config files will be stored here)

### Specify the following model options
# simulation_time            : simulation start/end times (example is given below)
# model_option               : model option (see available options below)
# surface_runoff_scheme      : surface runoff scheme for CFE and LASAM OPTION=[GIUH, NASH_CASCADE]
# precip_partitioning_scheme : preciptation partitioning scheme for CFE OPTION=[Schaake, Xinanjiang]
# is_netcdf_forcing          : True if forcing data is in netcdf format, otherwise False
# clean_all                  : delete all old/existing config files directories
# partition_pgkg             : Set to True for partitioning the geopackage for a parallel ngen run

#NOTE: The script assumes that .gpkg files are stored under cat_id/data/gage_id.gpkg (modify this part according to your settings)
#      the partitioned .json is also stored in the data directory
################################################ ###################################

root_dir     = "/Users/ahmadjan/Core/SimulationsData/worktasks/CAMELS_2024/basins516/"
workflow_dir = "/Users/ahmadjan/Core/SimulationsData/preprocessing/hydrofabric/smp_basin_workflow/basin_workflow/generate_files"
forcing_dir  = "/Users/ahmadjan/Core/SimulationsData/worktasks/CAMELS_2024/forcings/"
ngen_dir     = "/Users/ahmadjan/codes/ngen/ngen"

# simulation time format YYYYMMDDHHMM (YYYY, MM, DD, HH, MM)
simulation_time            = '{"start_time" : "2010-10-01 00:00:00", "end_time" : "2015-10-01 00:00:00"}' 
model_option               = "NC"
precip_partitioning_scheme = 'Schaake'
surface_runoff_scheme      = 'NASH_CASCADE'
is_netcdf_forcing          = True
clean_all                  = True
partition_gpkg             = True
print (simulation_time)

############ CHECKS ###################
assert (os.path.exists(root_dir))
assert (os.path.exists(workflow_dir))
assert (os.path.exists(forcing_dir))
assert (os.path.exists(ngen_dir))
######################################

gpkg_dirs = glob.glob(root_dir + "/*/", recursive = True)

if (is_netcdf_forcing):
    forcing_files = glob.glob(forcing_dir + "*.nc", recursive = True)
else:
    forcing_files = forcing_dir


#################################################################################
############################### MAIN LOOP #######################################
#################################################################################

for dir in gpkg_dirs:
    os.chdir(dir)
    
    print ("cwd: ", os.getcwd())
    gpkg_name     = os.path.basename(glob.glob(dir + "/data/*.gpkg")[0])  # <---- modify this line according to your settings
    gpkg_dir      = f"data/{gpkg_name}"                                   # <---- modify this line according to your settings
    
    print("=========================================")
    print ("Running : ", gpkg_dir)
    
    output_dir    = "inputs" + model_option

    if (os.path.isdir(output_dir) and clean_all):
        shutil.rmtree(output_dir)
        os.mkdir(output_dir)
    elif (not os.path.isdir(output_dir)):
        os.mkdir(output_dir)

    if (os.path.isdir("dem")):
        shutil.rmtree("dem")

    # find forcing data for the geopackage
    if(is_netcdf_forcing):
        id =  int(gpkg_name[:-5].split("_")[1]) # -5 is to remove .gpkg from the string
        forcing_file = [f for f in forcing_files if str(id) in f]
        if(len(forcing_file) == 1):
            forcing_dir = forcing_file[0]
        else:
            continue
    assert (os.path.exists(forcing_dir))
    
    workflow_driver = os.path.join(workflow_dir,"driver.py")

    driver = f'python {workflow_driver} -gpkg {gpkg_dir} -ngen {ngen_dir} -f {forcing_dir} \
    -o {output_dir} -m {model_option} -p {precip_partitioning_scheme} -r {surface_runoff_scheme} -t \'{simulation_time}\' \
    -netcdf {is_netcdf_forcing}'

    result = subprocess.call(driver,shell=True)


    #####################################################################
    # Parition geopackage
    if(partition_gpkg ):
        if (os.path.exists(f"{ngen_dir}/cmake_build/partitionGenerator")):
            
            x = gpd.read_file(gpkg_dir, layer="divides")
            num_div = len(x["divide_id"])
            nproc = 2
            if (num_div >= 4 and num_div <= 16):
                nproc = 4
            elif(num_div <= 48):
                nproc = 8
            elif(num_div <= 96):
                nproc = 16
            else:
                nproc = 20
            fpar = os.path.join("data", gpkg_name[:-5].split(".")[0] + f"-par{nproc}.json")     # -5 is to remove .gpkg from the string
            partition=f"{ngen_dir}/cmake_build/partitionGenerator {gpkg_dir} {gpkg_dir} {fpar} {nproc} \"\" \"\" "
            result = subprocess.call(partition,shell=True)
            
    #break
    
"""
coupled_models_options = {
"C"   : "cfe",
"L"   : "lasam",
"NC"  : "nom_cfe",
"NL"  : "nom_lasam",
"NCSS": "nom_cfe_smp_sft",
"NLSS": "nom_lasam_smp_sft",
"NT"  : "nom_topmodel",
"BC"  : "baseline_cfe",
"BL"  : "baseline_lasam"
}
"""
