## Installation instructions
Detailed instructions on how to build SoilMoistureProfiles for different setups/modes (standalone, with topmodel, and nextgen framework) are provided below. The schemes assume you have [GCC](https://gcc.gnu.org) and [CMAKE](https://cmake.org/) on your machine. These are simple examples providing a one-way coupling (model --> SoilMoistureProfiles) and are only for demonstration purposes. However, in real examples when coupled to other models, the soil moisture profiles effect freeze-thaw and other hydrological processes (such as infiltration).

### Build (standalone)
 ```
 git clone https://github.com/NOAA-OWP/SoilMoistureProfiles && cd SoilMoistureProfiles
 mkdir build && cd build
 cmake ../ -DSTANDALONE=ON
 make && cd ..
 ./run_smp.sh STANDALONE 
 ```

### Build (pseudo framework)
 ```
 git clone https://github.com/NOAA-OWP/SoilMoistureProfiles && cd SoilMoistureProfiles 
 git clone https://github.com/NOAA-OWP/topmodel extern/topmodel
 mkdir build && cd build
 cmake ../ -DWITHTOPMODEL=ON
 make && cd ..
 ./run_smp.sh WITHTOPMODEL] 
 ```
 
 ### Build (nextgen framework)
 Detailed instructions for running and building SoilMoistureProfiles coupled to other models (for instance, CFE or SFT) in the nextgen framework are provided at [instructions](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/INSTALL.md).
 