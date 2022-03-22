#!/bin/bash
g++ -lm -Wall -O -g ./src/bmi_main_mapper.cxx ./src/bmi_mapper.cxx ./src/smc_mapping.cxx -o run_bmi_mapper
#${CXX} -lm -Wall -O -g ./src/bmi_main.cxx ./src/bmi_coupler.cxx ./src/smc_profile.cxx -o run_bmi_mapper
