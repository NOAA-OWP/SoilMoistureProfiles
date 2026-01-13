#!/bin/bash
if [[ -n "${CXX:-}" ]] && command -v "$CXX" >/dev/null 2>&1; then
    : # CXX is set and valid, do nothing
elif command -v g++ >/dev/null 2>&1; then
    export CXX=g++
elif command -v clang++ >/dev/null 2>&1; then
    export CXX=clang++
else
    echo "ERROR: No C++ compiler found (CXX, g++, or clang++)" >&2
    exit 1
fi

"${CXX}" ./main_unittest.cxx ../src/bmi_soil_moisture_profile.cxx ../src/soil_moisture_profile.cxx -Wall -O -g -o run_smp -lm -lboost_serialization

if [[ $? -ne 0 ]]; then
    echo "ERROR: Compilation failed" >&2
    exit 1
fi

./run_smp configs/unittest_conceptual.txt configs/unittest_layered_constant.txt configs/unittest_layered_linear.txt
rm run_smp
rm -rf run_smp.dSYM
