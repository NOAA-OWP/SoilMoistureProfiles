cmake_minimum_required(VERSION 3.10)

set(CMAKE_C_COMPILER $ENV{CC})
set(CMAKE_CXX_COMPILER $ENV{CXX})


# set the project name
project(smp)

set(CMAKE_BUILD_TYPE Debug)
IF(CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug build.")
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)
message(CMAKE_CXX_COMPILER " ${CMAKE_CXX_COMPILER}")
message(CMAKE_C_COMPILER " ${CMAKE_C_COMPILER}")
message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")

# add the executable

#add_executable(smp ./src/bmi_main.cxx ./src/soil_moisture_profile.cxx ./src/bmi_coupler.cxx)

# ngen (add shared library
add_library(smpbmi SHARED ./src/bmi_coupler.cxx ./src/soil_moisture_profile.cxx ./include/bmi_coupler.hxx ./include/soil_moisture_profile.hxx)

#target_link_libraries(smp LINK_PUBLIC smplib)
#target_link_libraries(cfesft PRIVATE m)

target_include_directories(smpbmi PRIVATE /include)

##cmake -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_C_COMPILER=$CC -DCMAKE_INSTALL_PREFIX=`pwd` -DCMAKE_BUILD_TYPE=Debug ../
##cmake --build .



