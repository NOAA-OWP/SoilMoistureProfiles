#!/bin/bash
${CXX} -lm -Wall -O -g ./main_unit_test_bmi.cxx ../src/bmi_soil_moisture_profile.cxx ../src/soil_moisture_profile.cxx -o run_smp
./run_smp configs/unittest.txt
rm run_smp
rm -rf run_smp.dSYM
