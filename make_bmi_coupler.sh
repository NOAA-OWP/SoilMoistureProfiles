#!/bin/bash
g++ -lm -Wall -O -g ./src/bmi_main_moisture.cxx ./src/bmi_coupler.cxx ./src/soil_moisture_profile.cxx -o run_bmi_coupler
#${CXX} -lm -Wall -O -g ./src/bmi_main.cxx ./src/bmi_coupler.cxx ./src/smc_profile.cxx -o run_bmi_coupler
