# Installation instructions
Detailed instructions on how to build and run SoilMoistureProfiles for different setups/modes (pseudo framework with a conceptual soil reservoir, pseudo framework with soil moisture deficit, and nextgen framework) are provided below. The schemes assume you have [GCC](https://gcc.gnu.org) and [CMAKE](https://cmake.org/) on your machine. These are simple examples providing a one-way coupling (model --> SoilMoistureProfiles) and are only for demonstration purposes. However, in real examples when coupled to other models, the soil moisture profiles effect freeze-thaw and other hydrological processes (such as infiltration).

**Note**: Before running the following examples, it is recommended to run the unittests [tests](https://github.com/NOAA-OWP/SoilMoistureProfiles/tree/ajk/doc_update/tests).

## Pseudo framework example (CFE)
This example couples SMP to the [Conceptual Functional Equivalent (CFE)](https://github.com/NOAA-OWP/cfe) model which uses a conceptual soil reservoir.  SMP is used to compute `watertable` and `soil_moisture_profile`.
### Build
 ```
 git clone https://github.com/NOAA-OWP/SoilMoistureProfiles && cd SoilMoistureProfiles
 git clone https://github.com/NOAA-OWP/cfe extern/cfe
 mkdir build && cd build
 cmake ../ -DCFE=ON
 make && cd ..
 ```

### Run
<pre>
Run: <a href="https://github.com/NOAA-OWP/SoilMoistureProfiles/blob/ajk/doc_update/run_sft.sh">./run_smp.sh</a> CFE (from SoilMoistureProfiles directory)
</pre>

## Pseudo framework example (TopModel)
This example couples [TopModel](https://github.com/NOAA-OWP/topmodel), which tracks soil moisture deficit, to SMP.  SMP is used to compute `watertable` and the `soil_moisture_profile` from the soil moisture deficit.
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
The [nextgen framework](https://github.com/NOAA-OWP/ngen) allows the user to easily couple different BMI enabled models or model components, including SMP.  Detailed instructions for running and building SoilMoistureProfiles coupled to other models (for instance, CFE or [SoilFreezeThaw](https://github.com/NOAA-OWP/SoilFreezeThaw)) in the nextgen framework are provided at [instructions](https://github.com/NOAA-OWP/SoilFreezeThaw/blob/master/INSTALL.md).  Once nextgen is built, refer to the commands below to run the SMP examples in nextgen.
### Build
```
mkdir smp && cd smp (in the nextgen directory)
ln -s ../extern
ln -s ../data
ln -s extern/SoilMoistureProfiles/SoilMoistureProfiles/realizations
cp extern/topmodel/topmodel/data/* data
```
### Run
 - SMP with CFE example
   ```
   ../cmake_build/ngen data/catchment_data.geojson cat-27 data/nexus_data.geojson nex-26 realizations/realization_config_smp.json
   ```
 - SMP with Topmodel example
    ```
   cp ./extern/topmodel/topmodel/data/* data
   ../cmake_build/ngen data/catchment_data.geojson cat-27 data/nexus_data.geojson nex-26 realizations/realization_config_smp_topmodel.json
   ```
