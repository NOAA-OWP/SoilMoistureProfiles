############################################################################################
# Author  : Ahmad Jan
# Contact : ahmad.jan@noaa.gov
# Date    : October 11, 2023 
############################################################################################


import os,sys
import subprocess
import pandas as pd

# Note #1: from the command line just run 'python path_to/main.py'
# Note #2: make sure to adjust the following required arguments
# Note #3: several model coupling options are available, the script currently supports a few of them, for full list see
#          below the coupled_models_option


## Required arguments
# -gpkg : hydrofabric geopackage file (.gpkg)
# -f    : path to lumped forcings directory (pre-downloaed forcing data for each catchment)
# -ngen : path to nextgen directory (recommended to have a soft link in the current directory, ln -s path_to_ngen)
# -o    : path to output directory (all config files will be stored here)
# -m    : model option (see available options below)
# -p    : preciptation partitioning scheme for CFE OPTION=[Schaake, Xinanjiang]
# -r    : surface runoff scheme for CFE and LASAM OPTION=[GIUH, NASH_CASCADE]
# -t    : simulation start/end times (example is given below)

path_gpkg      = "data/Gage_05593900.gpkg"
path_forcing   = "data/forcing"
path_ngen      = "../ngen/"
path_output    = "inputsNC"
model_option   = "NC"
precip_partitioning_scheme  = 'Schaake'
surface_runoff_scheme  = 'NASH_CASCADE'



"""
# SCAN SITES

huc12_id = ['HUC12-100302030705', 'HUC12-101800090503', 'HUC12-100401041801', 'HUC12-102600120103',
            'HUC12-101301020105', 'HUC12-160101010403']
path_gpkg = f"data/{huc12_id[1]}.gpkg"
path_forcing = "/Users/ahmadjan/Core/simulations/owp/AGU_2023/forcing_all"
path_ngen = "../ngen_py3.11"
path_output = "inputsC"
model_option = "NCSS" #"BC"
#model_option = "NLSS"
runoff_scheme= 'Schaake'
#runoff_scheme= 'Xinanjiang'
"""


path_crf = os.path.dirname(sys.argv[0])

simulation_time = '{"start_time" : "2010-10-01 00:00:00", "end_time" : "2015-10-01 00:00:00"}' # format YYYYMMDDHHMM (YYYY, MM, DD, HH, MM)

print (simulation_time)

print ("CRF workflow path: ", path_crf )

path_crf_create = os.path.join(path_crf,"driver.py")

driver = f'python {path_crf_create} -gpkg {path_gpkg} -ngen {path_ngen} -f {path_forcing} \
-o {path_output} -m {model_option} -p {precip_partitioning_scheme} -r {surface_runoff_scheme} -t \'{simulation_time}\' '

print ("Running (from main.py):\n",driver)

result = subprocess.call(driver,shell=True)
    
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
