# Installation instructions
Detailed instructions on how to build and run SoilMoistureProfiles for different setups/modes (standalone, with topmodel, and nextgen framework) are provided below. The schemes assume you have [GCC](https://gcc.gnu.org) and [CMAKE](https://cmake.org/) on your machine. These are simple examples providing a one-way coupling (model --> SoilMoistureProfiles) and are only for demonstration purposes. However, in real examples when coupled to other models, the soil moisture profiles effect freeze-thaw and other hydrological processes (such as infiltration).

**Note**: Before running the following examples, it is recommended to run the unittests [tests](https://github.com/NOAA-OWP/SoilMoistureProfiles/tree/ajk/doc_update/tests).

## Pseudo framework example (CFE)
The example couples SMP to CFE (a conceptual soil reservoir) to compute `watertable` and `soil_moisture_profile`.
### Build
 ```
 git clone https://github.com/NOAA-OWP/SoilMoistureProfiles && cd SoilMoistureProfiles
 git clone https://github.com/NOAA-OWP/cfe extern/cfe
 mkdir build && cd build
 cmake ../ --DCFE=ON
 make && cd ..
 ```

### Run
<pre>
Run: <a href="https://github.com/NOAA-OWP/SoilMoistureProfiles/blob/ajk/doc_update/run_sft.sh">./run_smp.sh</a> CFE (from SoilMoistureProfiles directory)
</pre>

## Pseudo framework example (TopModel)
An example coupling TopModel to SMP to compute `watertable` and `soil_moisture_profile` from soil moisture deficit.
### Build
 ```
 git clone https://github.com/NOAA-OWP/SoilMoistureProfiles && cd SoilMoistureProfiles
 git clone https://github.com/NOAA-OWP/topmodel extern/topmodel
 mkdir build && cd build
 cmake ../ -DTOPMODEL=ON
 make && cd ..
 ```
### Run
<pre>
Run: <a href="https://github.com/NOAA-OWP/SoilMoistureProfiles/blob/ajk/doc_update/run_sft.sh"> ./run_smp.sh</a> TOPMODEL (from SoilMoistureProfiles directory)
</pre>

## Nextgen framework example
Detailed instructions for running and building SoilMoistureProfiles coupled to other models (for instance, CFE or SFT) in the nextgen framework are provided at [instructions](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/INSTALL.md).
### Build
```
mkdir smp && cd smp (in the nextgen directory)
ln -s ../extern
ln -s ../data
ln -s extern/SoilMoistureProfiles/SoilMoistureProfiles/realizations
cp extern/topmodel/topmodel/data/* data
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
