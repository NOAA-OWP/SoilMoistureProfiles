cmake_minimum_required(VERSION 3.12)

# Uncomment this and rebuild artifacts to enable debugging
set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(NGEN "NGEN" OFF)
option(MODELC "MODELC" OFF)
option(MODELCXX "MODELCXX" OFF)

message(CMAKE_CXX_COMPILER " ${CMAKE_CXX_COMPILER}")
message(CMAKE_C_COMPILER " ${CMAKE_C_COMPILER}")
message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

if(NGEN)
  message("ngen framework build!")
  add_definitions(-DNGEN)
endif()

if(MODELC)
 set(exe_name "model_c")
 add_definitions(-DMODELC)
 message("Model C build!")
endif()

if(MODELCXX)
 set(exe_name "model_cxx")
 add_definitions(-DMODELCXX)
 message("Model CXX build!")
endif()

#set the project name
project(bmi_models VERSION 1.0.0 DESCRIPTION "OWP BMI wrappers for C and CXX models")


if(MODELC)
  add_executable(${exe_name} ./src/main.cxx ./src/bmi_model_c.cxx ./src/model_c.cxx)
  target_link_libraries(${exe_name} PRIVATE m)
elseif(MODELCXX)
  add_executable(${exe_name} ./src/main.cxx ./src/bmi_model_cxx.cxx ./src/model_cxx.cxx)
  target_link_libraries(${exe_name} PRIVATE m)
endif()


unset(MODELC CACHE)
unset(MODELCXX CACHE)
unset(NGEN CACHE)
