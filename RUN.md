# Running instructions
Here, we provide a few examples to run and test SMP (Soil Moisture Profiles). Before running the following examples, it is recommended to run the unittest [tests](https://github.com/NOAA-OWP/SoilMoistureProfiles/tree/ajk/doc_update/tests).

### Standalone example
  - the example runs SMP for a hypothetical conceptual soil reservoir (with prescribed soil_storage and soil_storage_change SMP BMI input variables) to compute `watertable` and `soil_moisture_profile` and compare the results against benchmark results. The hypothetical conceptual soil reservoir mimics CFE.
<pre>
Run: <a href="https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/doc_update/run_sft.sh">./run_smp.sh</a> STANDALONE (from SoilMoistureProfiles directory)    
</pre>

### Pseudo framework example
 - An example coupling TopModel to SMP to compute `watertable` and `soil_moisture_profile` from soil moisture deficit.

<pre>
Run: <a href="https://github.com/NOAA-OWP/SoilFreezeThaw/blob/ajk/doc_update/run_sft.sh"> ./run_smp.sh</a> WITHTOPMODEL (from SoilMoistureProfiles directory)    
</pre>

