# Running instructions
Here, we provide a few examples to run and test SMP (Soil Moisture Profiles). Before running the following examples, it is recommended to run the unittest [tests](https://github.com/NOAA-OWP/SoilMoistureProfiles/tree/ajk/doc_update/tests).

### Standalone example
  - the example runs SMP for a hypothetical conceptual soil reservoir (with prescribed soil_storage and soil_storage_change SMP BMI input variables) to compute `watertable` and `soil_moisture_profile` and compare the results against benchmark results. The hypothetical conceptual soil reservoir mimics CFE.
<pre>
Run: <a href="https://github.com/NOAA-OWP/SoilMoistureProfiles/blob/ajk/doc_update/run_sft.sh">./run_smp.sh</a> STANDALONE (from SoilMoistureProfiles directory)    
</pre>

### Pseudo framework example
 - An example coupling TopModel to SMP to compute `watertable` and `soil_moisture_profile` from soil moisture deficit.

<pre>
Run: <a href="https://github.com/NOAA-OWP/SoilMoistureProfiles/blob/ajk/doc_update/run_sft.sh"> ./run_smp.sh</a> WITHTOPMODEL (from SoilMoistureProfiles directory)    
</pre>

### Nextgen framework example
```
mkdir smp && cd smp (in the nextgen directory)
ln -s ../extern
ln -s ../data
ln -s extern/SoilMoistureProfiles/SoilMoistureProfiles/realizations

 - Standalone example
   ../cmake_build/ngen data/catchment_data.geojson cat-27 data/nexus_data.geojson nex-26 realizations/realization_config_smp.json 
 - With Topmodel example
   cp ./extern/topmodel/topmodel/data/* data
   ../cmake_build/ngen data/catchment_data.geojson cat-27 data/nexus_data.geojson nex-26 realizations/realization_config_smp_topmodel.json (for with topmodel example)
```

