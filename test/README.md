# Unit test for Soil Moisture Profiles model
Usage: run `./run_unittest.sh` (replace/set $CXX to your g++ compiler)

 - Multiple checks are performed to test the functionality of the BMI functions:
   - Check number of input/output variables
   - Test `GetVar*` methods for the input variables and compare against initial data; if failed, will throw an error
   - Test `GetGrid*` methods and compare data against known values.
   - Loop over the input variables, use `Set*` and `Get*` methods to verify `Get*` return the same data set by `Set*`
   - Step (4) for output variables
   - Using `Update` method, advance the model to get updated soil moisture profile and water table location. Compare the `soil moisture profile` and `water table location` against benchmark values
 - Conceptual Soil Reservoir unittest
   - the test uses an input `soil_storage` (to mimic CFE input) to compute soil moisture profiles and compare against benchmark
 - Layered Soil Reservoir unittest
   - Test scenario #1: four wetting fronts with lower two layers fully saturated; and two wetting fronts in the 2nd layer (CONSTANT PROFILE OPTION)
   - Test scenario #2: four wetting fronts with only top wetting fully saturated; two wetting fronts in the top layer (CONSTANT PROFILE OPTION)
   - Test scenario #3: four wetting fronts with only top wetting fully saturated; two wetting fronts in the top layer (LINEAR PROFILE OPTION)