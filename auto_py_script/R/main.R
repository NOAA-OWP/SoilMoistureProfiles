# @author Ahmad Jan Khattak
# @email ahmad.jan@noaa.gov
# @date  December 22, 2023

# The script downloads a geopackge given USGS gauge id, computes TWI and GIUH
# Appends model attributes from the hydrofabric parquet file to the geopackage

######################## REQUIRED INPUT #######################################
# The script downloads geopackage given a gage ID or the user provides the basin 
# gpkg file  (SEE STEP #3). The script is divided into few steps:
# STEP #1: Setup/install packages needed for running the script
# STEP #2: Set pointers to the right directory (e.g., r_path, gpkg_path; see below)
# STEP #3: Download geopackage
# STEP #4: Add model attributes to the geopackage
# STEP #5: Computes and appends TWI and width function to the geopackage
# STEP #6: Computes and appends GIUH the geopackage
###############################################################################


######################### INSTALL REQUIRED PACKAGES ############################
# STEP #1: The packages need to run the script
################################################################################
if(!requireNamespace("hydrofabric", quietly=TRUE)) 
  devtools::install_github("noaa-owp/hydrofabric", ref = 'b07c109', force = TRUE)

if(!requireNamespace("climateR", quietly=TRUE)) 
  devtools::install_github("mikejohnson51/climateR", force = TRUE)

if(!requireNamespace("zonal", quietly=TRUE))
  devtools::install_github("mikejohnson51/zonal", ref = '24c5ee6', force = TRUE)

if(!requireNamespace("AOI", quietly=TRUE))
  devtools::install_github("mikejohnson51/AOI")


if(!requireNamespace("sf", quietly=TRUE)) install.packages("sf")
if(!requireNamespace("terra", quietly=TRUE)) install.packages("terra")
if(!requireNamespace("whitebox", quietly=TRUE)) install.packages("whitebox")
if(!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr")
if(!requireNamespace("glue", quietly=TRUE)) install.packages("glue")
if(!requireNamespace("raster", quietly=TRUE)) install.packages("raster")
if(!requireNamespace("jsonlite", quietly=TRUE)) install.packages("jsonlite")
if(!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")


library(hydrofabric)
library(climateR)
library(zonal)
library(whitebox)
library(sf)
library(terra)
library(dplyr)
library(glue)
library(raster)
library(jsonlite)
library(ggplot2)


############################## SET PATHS ######################################
# STEP #2. Note the three substeps (a), (b), and (c)
################################################################################

# (a) Point r_path to the directory of R scripts downloaded from the repository
r_path = "/Users/ahmadjan/Core/SimulationsData/preprocessing/hydrofabric/R"
source(glue("{r_path}/twi_width_function.R"))
source(glue("{r_path}/helper.R"))
source(glue("{r_path}/giuh_function.R"))

# (b) Point gpkg_path to the directory where geopackage and other related files will be stored
# NOTE: .gpkg file should be stored under {gpkg_path}/data directory
gpkg_path = "/Users/ahmadjan/Core/SimulationsData/preprocessing/hydrofabric/sebecx"
setwd(gpkg_path)
wbt_wd(getwd())

# (c) DEM (located at S3 endpoint)
dem_path = "/vsicurl/https://lynker-spatial.s3.amazonaws.com/gridded-resources/dem.vrt"

# DEM and related files (such as projected/corrected DEMs, and specific contributing area rasters are stored here)
directory="dem"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)




######################### DOWNLOAD GEOPACKAGE ##################################
# STEP #3: provide USGS gauge id or your own geopackage (set path in the ELSE block)
################################################################################
is_gpkg_provided = FALSE
if(!is_gpkg_provided) {
  gage_id = '01033000' # A hydrolocation URI
  hl = glue('Gages-{gage_id}') # sebec river basin in Maine
  outfile <- glue('data/gage_{gage_id}.gpkg')

  ## caching the downloaded VPU files to "data" and writing all layers to "outfile"
  subset_network(hl_uri = hl, cache_dir = "data", outfile = outfile)
} else {
  outfile <- glue('data/gage_{gage_id}.gpkg') # <------ set this path if .gpkg is store somewhere else
}

## Stop if .gpkg does not exist
if (!file.exists(outfile))
  stop(glue("FILE '{outfile}' DOES NOT EXIST!!"))

div <- read_sf(outfile, 'divides')
nexus <- read_sf(outfile, 'nexus')
streams <- read_sf(outfile, 'flowpaths')




########################## MODELS' ATTRIBUTES ##################################
# STEP #4: Add models' attributes from the parquet file to the geopackage
################################################################################
# print layers before appending model attributes
layers_before_cfe_attr <- sf::st_layers(outfile)
print (layers_before_cfe_attr$name)

m_attr <- add_model_attributes(div_path = outfile)

layers_after_cfe_attr <- sf::st_layers(outfile)
print (layers_after_cfe_attr$name)


############################### GENERATE TWI ##################################
# STEP #5: Generate TWI and width function and write to the geopackage
###############################################################################
dem_function(infile = outfile, directory, dem_path)

twi <- twi_function(infile = outfile, directory = directory, distribution = 'simple', nclasses = 30)

width_dist <- width_function(outfile, directory = directory)

dat_values = data.frame(ID = twi$divide_id, twi = twi$fun.twi, width_dist = width_dist$fun.downslope_fp_length)

# write TWI and width function layers to the geopackage
names(dat_values)
colnames(dat_values) <- c('divide_id', 'twi', 'width_dist')
names(dat_values)
sf::st_write(dat_values, outfile, layer = "twi", append = FALSE)

width_dist
############################### GENERATE GIUH ##################################
# STEP #6: Generate GIUH and write to the geopackage
################################################################################
# There are many "model" options to specify the velocity.
# Here we are using a simple approach: constant velocity as a function of upstream drainage area.
vel_channel     <- 1.0
vel_overland    <- 0.5
vel_gully       <- 0.2
gully_threshold <- 90.0

giuh_compute <- giuh_function(infile = outfile, directory, vel_channel,
                              vel_overland, vel_gully, gully_threshold)


# write GIUH layer to the geopackage
giuh_dat_values = data.frame(ID = giuh_compute$divide_id, giuh = giuh_compute$fun.giuh_minute)
names(giuh_dat_values)
colnames(giuh_dat_values) <- c('divide_id', 'giuh')
names(giuh_dat_values)
sf::st_write(giuh_dat_values, outfile, layer = "giuh", append = FALSE)

###############################################################################
