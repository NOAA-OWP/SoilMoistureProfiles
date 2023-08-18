## Installation instructions
Detailed instructions on how to build SoilMoistureProfiles for different setups/modes (standalone, with topmodel, and nextgen framework) are provided below. The schemes assume you have [GCC](https://gcc.gnu.org) and [CMAKE](https://cmake.org/) on your machine. These are simple examples providing a one-way coupling (model --> SoilMoistureProfiles) and are only for demonstration purposes. However, in real examples when coupled to other models, the soil moisture profiles effect freeze-thaw and other hydrological processes (such as infiltration).
 
 ```
 git clone https://github.com/NOAA-OWP/SoilMoistureProfiles.git
 cd SoilMoistureProfiles 
 git clone https://github.com/NOAA-OWP/topmodel extern/topmodel (needed for running with Topmodel, i.e., -DWITHTOPMODEL=ON)
 mkdir build && cd build
 cmake ../ [-DSTANDALONE=ON,-DWITHTOPMODEL=ON] (pick one option, e.g. `cmake ../ -DSTANDALONE=ON`)
 make && cd ..
 ./run_smp.sh [STANDALONE, WITHTOPMODEL] (pick one option) 
 ```

 ### Coupled mode (ngen and pseudo frameworks):
  * Coupling SoilMoistureProfiles to any module (for instance, CFE or SFT) **must** follow these [instructions](https://github.com/NOAA-OWP/SoilFreezeThaw) for building and running. **Note separate instructions are provided for building/running in the ngen framework on the  SoilFreezeThaw repo.**
  * Follow this example: [couple SMP with SFT](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/src/main_cfe_aorc_pet_ftm.cxx)