cmake_minimum_required(VERSION 3.12)

# Uncomment this and rebuild artifacts to enable debugging
set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(smpbmi VERSION 1.0.0 DESCRIPTION "OWP SMP BMI Module Shared Library")

option(NGEN "NGEN" OFF)

if(NGEN)
   add_definitions(-DNGEN)
endif()

# SMP
set(SMP_LIB_NAME_CMAKE smpbmi)
set(SMP_LIB_DESC_CMAKE "OWP SMP BMI Module Shared Library")

# Make sure these are compiled with this directive
add_compile_definitions(BMI_ACTIVE)

if(WIN32)
    add_library(smpbmi src/bmi_soil_moisture_profile.cxx src/soil_moisture_profile.cxx)
else()
    add_library(smpbmi SHARED src/bmi_soil_moisture_profile.cxx src/soil_moisture_profile.cxx)
endif()

target_include_directories(smpbmi PRIVATE include)

set_target_properties(smpbmi PROPERTIES VERSION ${PROJECT_VERSION})

set_target_properties(smpbmi PROPERTIES PUBLIC_HEADER ./include/bmi_soil_moisture_profile.hxx)

# Code requires minimum of C99 standard to compile
set_target_properties(smpbmi PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON)

include(GNUInstallDirs)

install(TARGETS smpbmi
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

configure_file(smpbmi.pc.in smpbmi.pc @ONLY)

install(FILES ${CMAKE_BINARY_DIR}/smpbmi.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
