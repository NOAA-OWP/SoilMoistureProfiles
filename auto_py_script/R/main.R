# @author Ahmad Jan Khattak
# @email ahmad.jan@noaa.gov
# @date  December 22, 2023

# The script downloads a geopackge given USGS gauge id, computes TWI and GIUH
# appends model attributes from the hydrofabric parquet file to the geopackage

devtools::install_github("mikejohnson51/climateR", force = TRUE)
devtools::install_github("noaa-owp/hydrofabric")#,  force = TRUE)
devtools::install_github("mikejohnson51/zonal", force = TRUE)
devtools::install_github("mikejohnson51/AOI")

library(sf)
library(terra)
library(whitebox)
library(zonal)
library(dplyr)
library(glue)
library(raster)
library(hydrofabric)
library(zonal)
library(dplyr)
library(jsonlite)
library(tmap)
tmap_mode("view")
library(zoom)

r_path = "/Users/ahmadjan/Core/SimulationsData/preprocessing/hydrofabric/R"
source(glue("{r_path}/twi_width_function.R"))
source(glue("{r_path}/add_model_attributes.R"))

gpkg_path = "/Users/ahmadjan/Core/SimulationsData/preprocessing/hydrofabric/sebecx"
setwd(gpkg_path)
wbt_wd(getwd())

dem_path = "/vsicurl/https://lynker-spatial.s3.amazonaws.com/gridded-resources/dem.vrt"


################################################################
name="dem"
working_dir="test"
directory=glue("{working_dir}/{name}")
print (directory)
dir.create(directory, recursive = TRUE, showWarnings = FALSE)

###########################################################
# Step 1: Download geopackage
###########################################################
gage_id = '01033000' # A hydrolocation URI
hl = glue('Gages-{gage_id}') # sebec river basin
outfile = glue('data/gage_{gage_id}.gpkg')

## caching the downloaded VPU files to "data" and writing all layers to "outfile"
subset_network(hl_uri = hl, cache_dir = "data", outfile = outfile)

div <- read_sf(outfile, 'divides')
nexus <- read_sf(outfile, 'nexus')
streams <- read_sf(outfile, 'flowpaths')



###########################################################
# Step 2: Add model attributes to the geopackage
###########################################################
# print layers before appending model attributes
layers_before_cfe_attr <- sf::st_layers(outfile)
print (layers_before_cfe_attr$name)

m_attr <- add_model_attributes(div_path = outfile)

layers_after_cfe_attr <- sf::st_layers(outfile)
print (layers_after_cfe_attr$name)



###########################################################
# Step3: Generate TWI and width function and write to the geopackage
###########################################################

twi <- twi_function(outfile, directory = directory, dem_path = dem_path, 
                          nclasses = 5)

width_dist <- width_function(outfile, directory = directory)

dat_values = data.frame(ID = twi$divide_id, twi = twi$fun.twi, width_dist = width_dist$fun.down_fp_length)

names(dat_values)
colnames(dat_values) <- c('divide_id', 'twi', 'width_dist')
names(dat_values)
sf::st_write(dat_values, outfile, layer = "twi", append = FALSE)

#r_list <- fromJSON(width_dist['fun.down_fp_length']$fun.down_fp_length[1])
#print (r_list)
#print (sum(r_list$frequency))


###########################################################
# Step4: Generate GIUH and write to the geopackage
###########################################################
