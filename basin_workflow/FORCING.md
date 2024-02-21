# Forcing data via NextGen Forcings Engine
	
The NextGen Forcing Engine source repository is located [here](https://github.com/jduckerOWP/ngen-forcing). 
Make sure to go through the README.md to familiarize yourself with the forcing engine. At the time of writing this note, 
the user needs to be on VPN (GFE) to download AORC data directly and run the following script. 
For those who don’t have direct access to the data will need to request Jason Ducker (the main developer/maintainer of the 
NextGen [Forcings Engine](https://github.com/NOAA-OWP/ngen-forcing]) to get the processed data for their region of interest.

## Data download instructions:

### Clone ngen-forcing and create conda environment
- git clone https://github.com/jduckerOWP/ngen-forcing ngen/extern/
- conda create -n ngen_forcing
- conda activate ngen_forcing
- conda install -c conda-forge  pybind11 xarray dask netCDF4 numpy gdal  geopandas pandas wget scipy shapely multiprocess
- pip install wget

### Build Exactextract and pybindings
- git clone https://github.com/jdalrym2/exactextract
- cd exactextract
- git checkout coverage-fraction-pybindings
- Optional:
  - mkdir cmake-build-release
  - cd cmake-build-release
  - cmake -DCMAKE_BUILD_TYPE=Release .. -DBUILD_DOC=OFF
  - make
  - sudo make install
  - ./catch_tests
-	python3 setup.py bdist
-	python3 setup.py install
-	python
  - Run the following to test if all libraries will load without any issues
    ```
    from exactextract import GDALDatasetWrapper, GDALRasterWrapper, Operation, MapWriter, FeatureSequentialProcessor, GDALWriter
    ``` 
### Run NextGen Forcings driver to download the data
- create data/forcing (inside your basin directory)
- python ngen/extern/ngen-forcing/NextGen_Lumped_Forcings_Driver/Run_NextGen_lumped_driver.py
- Note: Go through the README.md [here](https://github.com/jduckerOWP/ngen-forcing/blob/master/NextGen_Lumped_Forcings_Driver/README.md)
  to understand and adjust different directories according to your settings.
