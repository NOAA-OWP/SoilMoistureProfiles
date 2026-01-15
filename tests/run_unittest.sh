#!/bin/bash
######################################################################
# Script to compile and run the unit tests for the Soil Moisture Profile
# models.
## Usage: ./run_unittest.sh
######################################################################

# Capture the directory of this script so we can use relative paths
THIS_DIRECTORY="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Make sure we have a C++ compiler
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

# Compile the unit test executable. Link to the math library and Boost
# serialization library.
"${CXX}" "${THIS_DIRECTORY}/main_unittest.cxx" \
        "${THIS_DIRECTORY}/../src/bmi_soil_moisture_profile.cxx" \
        "${THIS_DIRECTORY}/../src/soil_moisture_profile.cxx" \
        -Wall -O -g -o "${THIS_DIRECTORY}/run_smp" \
        -lm \
        -lboost_serialization

if [[ $? -ne 0 ]]; then
    echo "ERROR: Compilation failed" >&2
    exit 1
fi

# Run the unit tests with different configurations
"${THIS_DIRECTORY}/run_smp" "${THIS_DIRECTORY}/configs/unittest_conceptual.txt" \
    "${THIS_DIRECTORY}/configs/unittest_layered_constant.txt" \
    "${THIS_DIRECTORY}/configs/unittest_layered_linear.txt"

# Clean up the executable and debug files
rm "${THIS_DIRECTORY}/run_smp"
rm -rf "${THIS_DIRECTORY}/run_smp.dSYM"
