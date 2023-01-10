#!/bin/bash
g++ -lm -Wall -O -g ./src/bmi_main_regridder.cxx ./src/bmi_regridder.cxx ./src/soil_moisture_regridding.cxx -o run_bmi_mapper
#${CXX} -lm -Wall -O -g ./src/bmi_main.cxx ./src/bmi_coupler.cxx ./src/smc_profile.cxx -o run_bmi_mapper
