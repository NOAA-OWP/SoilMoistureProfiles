# Installation instructions
Detailed instructions on how to build and run SoilMoistureProfiles for different setups/modes (standalone, with topmodel, and nextgen framework) are provided below. The schemes assume you have [GCC](https://gcc.gnu.org) and [CMAKE](https://cmake.org/) on your machine. These are simple examples providing a one-way coupling (model --> SoilMoistureProfiles) and are only for demonstration purposes. However, in real examples when coupled to other models, the soil moisture profiles effect freeze-thaw and other hydrological processes (such as infiltration).

**Note**: Before running the following examples, it is recommended to run the unittests [tests](https://github.com/NOAA-OWP/SoilMoistureProfiles/tree/ajk/doc_update/tests).

## Standalone example
The example runs SMP for a hypothetical conceptual soil reservoir (with prescribed soil_storage and soil_storage_change SMP BMI input variables) to compute `watertable` and `soil_moisture_profile` and compare the results against benchmark results. The hypothetical conceptual soil reservoir mimics CFE.
### Build
 ```
 git clone https://github.com/NOAA-OWP/SoilMoistureProfiles && cd SoilMoistureProfiles
 mkdir build && cd build
 cmake ../ -DSTANDALONE=ON
 make && cd ..
 ```

### Run
<pre>
Run: <a href="https://github.com/NOAA-OWP/SoilMoistureProfiles/blob/ajk/doc_update/run_sft.sh">./run_smp.sh</a> STANDALONE (from SoilMoistureProfiles directory)    
</pre>

## Pseudo framework example
An example coupling TopModel to SMP to compute `watertable` and `soil_moisture_profile` from soil moisture deficit.
### Build
 ```
 git clone https://github.com/NOAA-OWP/SoilMoistureProfiles && cd SoilMoistureProfiles 
 git clone https://github.com/NOAA-OWP/topmodel extern/topmodel
 mkdir build && cd build
 cmake ../ -DWITHTOPMODEL=ON
 make && cd ..
 ```
### Run
<pre>
Run: <a href="https://github.com/NOAA-OWP/SoilMoistureProfiles/blob/ajk/doc_update/run_sft.sh"> ./run_smp.sh</a> WITHTOPMODEL (from SoilMoistureProfiles directory)    
</pre>

## Nextgen framework example
Detailed instructions for running and building SoilMoistureProfiles coupled to other models (for instance, CFE or SFT) in the nextgen framework are provided at [instructions](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/INSTALL.md).
### Build
```
mkdir smp && cd smp (in the nextgen directory)
ln -s ../extern
ln -s ../data
ln -s extern/SoilMoistureProfiles/SoilMoistureProfiles/realizations
```
### Run
 - Standalone example
   ```
   ../cmake_build/ngen data/catchment_data.geojson cat-27 data/nexus_data.geojson nex-26 realizations/realization_config_smp.json
   ```
 - With Topmodel example
    ```
   cp ./extern/topmodel/topmodel/data/* data
   ../cmake_build/ngen data/catchment_data.geojson cat-27 data/nexus_data.geojson nex-26 realizations/realization_config_smp_topmodel.json
   ``` 
