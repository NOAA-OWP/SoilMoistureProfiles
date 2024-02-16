############################################################################################
# Author  : Ahmad Jan
# Contact : ahmad.jan@noaa.gov
# Date    : October 11, 2023 
############################################################################################

# driver of the script, generates configuration files and realization files by calling
# the respective scripts

import os
import sys
import subprocess
from pathlib import Path
from typing import List, Union
import argparse
import json

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

class colors:
    BLUE    = '\33[34m'
    BLACK   = '\33[30m'
    RED     = '\33[31m'
    CYAN    = '\033[96m'
    GREEN   = '\033[32m'
    WARNING = '\033[93m'
    FAIL    = '\033[91m'
    ENDC    = '\033[0m'
    BOLD    = '\033[1m'
    UNDERLINE = '\033[4m'
    
def main():

    try:
        parser = argparse.ArgumentParser()
        parser.add_argument("-gpkg", dest="gpkg_file",     type=str, required=True,  help="the gpkg file")
        parser.add_argument("-f",    dest="forcing_dir",   type=str, required=True, help="the forcing files directory")
        parser.add_argument("-o",    dest="output_dir",    type=str, required=True,  help="the output files directory")
        parser.add_argument("-ngen", dest="ngen_dir",      type=str, required=True,  help="the ngen directory")
        parser.add_argument("-r",    dest="runoff_scheme", type=str, required=True, help="option for runoff scheme")
        parser.add_argument("-t",    dest="time",          type=str, required=True,  help="simulation start/end time")
        parser.add_argument("-m",    dest="models_option", type=str, required=True,  help="option for models coupling\n%s"%coupled_models_options)

        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    # check if geopackage file exists, if not throw an error
    if (not os.path.exists(args.gpkg_file)):
        sys.exit('The gpkg file (%s) does not exist!'%args.gpkg_file)

    # check if forcing dir exists, if not throw an error
    if (not os.path.exists(args.forcing_dir)):
        sys.exit('The forcing directory (%s) does not exist!'%args.forcing_dir)

    # check if output dir (to which output will be written to) exists, if not throw an error, user must create one 
    if (not os.path.exists(args.output_dir)):
        sys.exit("output dir does not exist, create one!")

    # not all model coupling options are support yet, throw an error if un-supported options are provided
    if (args.models_option in ["C","L","B"]):
        str_model = "Option under development: "+ str(coupled_models_options[args.models_option])+ ", quitting..."
        sys.exit(str_model)

    # check if the model option provided is valid and supported
    try:
        coupled_models = coupled_models_options[args.models_option]
    except:
        parser.print_help()
        str_msg = "*** Invalid model option provided: (%s) ***"%args.models_option
        sys.exit(str_msg)

    path_crf = os.path.dirname(sys.argv[0])

    # Note: for baseline simulations, models coupling is still either NCSS or NLSS,
    #       only realization file changes to add jinjaBMI, more output vars etc.
    baseline_case = False
    if ( 'baseline_cfe' == coupled_models_options[args.models_option]):
        baseline_case = True
        coupled_models = coupled_models_options["NCSS"]
    elif ( 'baseline_lasam' == coupled_models_options[args.models_option]):
        baseline_case = True
        coupled_models = coupled_models_options["NLSS"]


    path_crf_gen_files = os.path.join(path_crf,"configuration.py")

    generate_config_files = f'python {path_crf_gen_files} -gpkg {args.gpkg_file} -ngen {args.ngen_dir} \
                              -f {args.forcing_dir} -o {args.output_dir} -m {coupled_models} -r {args.runoff_scheme} -t \'{args.time}\' '

    print ("*******************************************")
    print (colors.BLUE)
    print ("Running (from driver.py):\n", generate_config_files)
    print ("Model option provided: ", args.models_option)
    print ("Generating configuration files for model(s) option: ", coupled_models_options[args.models_option])
    print (colors.ENDC)
    print ("*******************************************")
    

    # the following call will generate config files
    result = subprocess.call(generate_config_files,shell=True)

    # if the results is FALSE, something went wrong while generating config files, so throw an error
    if (result):
        sys.exit("config files could not be generated, check the options provided!")

    print ("*******************************************")
    print (colors.GREEN)
    print ("Generating realization file ...")

    path_crf_real_file = os.path.join(path_crf,"realization.py")

    generate_realization_file = f'python {path_crf_real_file} -ngen {args.ngen_dir} -f {args.forcing_dir} -i {args.output_dir} \
                                  -m {coupled_models} -r {args.runoff_scheme} -b {baseline_case} -t \'{args.time}\' '

    print ("Running (from driver.py): \n ", generate_realization_file)

    print (colors.ENDC)
    
    result = subprocess.call(generate_realization_file,shell=True)

    
    if (result):
        sys.exit("realization file could not be generated, check the options provided!")
    else:
        print ("************* DONE ************** ")

    if (baseline_case):
        infile = os.path.join(os.getcwd(), "realization_%s.json"%coupled_models)
        outfile = os.path.join(os.getcwd(), "realization_%s_baseline.json"%coupled_models)
        
        path_crf_baseline_file = os.path.join(path_crf,"baseline.py")
        
        convert_to_baseline = f'python {path_crf_baseline_file} -i {infile} -idir {args.output_dir} -o {outfile} -m {coupled_models}'

        print (colors.BLUE)
        print ("Generating baseline realization file ...")
        print ("Running (from driver.py):\n ", convert_to_baseline)
        print (colors.ENDC)
        
        result = subprocess.call(convert_to_baseline,shell=True)

        if (result):
            sys.exit("realization file could not be generated, check the options provided!")
        else:
            print ("************* DONE (Baseline realization file successfully generated!) ************** ")
        
if __name__ == "__main__":
    main()

