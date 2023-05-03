#!/bin/bash
${CXX} -lm -Wall -O -g ./main_unittest.cxx ../src/bmi_soil_moisture_profile.cxx ../src/soil_moisture_profile.cxx -o run_smp
./run_smp configs/unittest_conceptual.txt configs/unittest_layered_constant.txt configs/unittest_layered_linear.txt
rm run_smp
rm -rf run_smp.dSYM
