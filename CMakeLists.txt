cmake_minimum_required(VERSION 3.12)

string(ASCII 27 Esc)
set(ColourReset "${Esc}[m")
set(Red         "${Esc}[31m")

# Uncomment this and rebuild artifacts to enable debugging
set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")

project(smpbmi VERSION 1.0.0 DESCRIPTION "OWP SMP BMI Module Shared Library")

option(NGEN "NGEN" OFF)
option(CFE "CFE" OFF)
option(TOPMODEL "TOPMODEL" OFF)

if(NGEN)
   add_definitions(-DNGEN)
endif()

if(CFE)
message("${Red} Soil moisture profile for CFE! ${ColourReset}")
set(exe_name "smp_cfe")
elseif(TOPMODEL)
message("${Red} Soil moisture profile for topmodel (watertable-based)! ${ColourReset}")
set(exe_name "smp_topmodel")
endif()

# SMP
set(SMP_LIB_NAME_CMAKE smpbmi)
set(SMP_LIB_DESC_CMAKE "OWP SMP BMI Module Shared Library")


# Make sure these are compiled with this directive
add_compile_definitions(BMI_ACTIVE)

# add the executable
if(CFE)
  add_executable(${exe_name} ./src/main_smp_cfe.cxx ./extern/cfe/src/cfe.c ./extern/cfe/src/bmi_cfe.c
                             ./extern/cfe/src/giuh.c ./extern/cfe/src/conceptual_reservoir.c ./extern/cfe/src/nash_cascade.c)
 target_include_directories(${exe_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/extern/cfe/include)
elseif(TOPMODEL)
  add_executable(${exe_name} ./src/main_smp_topmodel.cxx ./extern/topmodel/src/topmodel.c ./extern/topmodel/src/bmi_topmodel.c)
  target_include_directories(${exe_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/extern/topmodel/include)
endif()


if(TOPMODEL OR CFE)
  add_library(smplib ./src/bmi_soil_moisture_profile.cxx ./src/soil_moisture_profile.cxx ./include/bmi_soil_moisture_profile.hxx
                     ./include/soil_moisture_profile.hxx)
  target_link_libraries(${exe_name} LINK_PUBLIC smplib)
  target_include_directories(${exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
  target_link_libraries(${exe_name} PRIVATE m)
endif()


## for NGEN BUILD
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

unset(STANDALONE CACHE)
unset(WITHTOPMODEL CACHE)
