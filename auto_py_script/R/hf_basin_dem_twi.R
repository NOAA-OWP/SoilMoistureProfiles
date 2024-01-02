# install ClimateR if running this for the first time
devtools::install_github("mikejohnson51/climateR", force = TRUE)
devtools::install_github("noaa-owp/hydrofabric")#,  force = TRUE)
devtools::install_github("mikejohnson51/zonal", force = TRUE)
devtools::install_github("mikejohnson51/AOI")
ARROW_R_DEV=TRUE
#install.packages("arrow")

library(dplyr)
library(arrow)
library(hydrofabric)
library(zonal)
library(raster)
library(jsonlite)
#library(climateR);
library(sf); 
library(terra)
library(mapview)
library(AOI)
##########################################################################
# set root path
r_path = "/Users/ahmadjan/Core/SimulationsData/preprocessing/hydrofabric/sebecx"
setwd(r_path)

# A hydrolocation URI
gage_id = '01033000'
hl = glue('Gages-{gage_id}') # sebec river basin
outfile = glue('data/gage_{gage_id}.gpkg')
outfile_bbox = glue('data/gage_{g_id}-bbox.gpkg')
# The output directory (the geopackage will be stored in data directory)
print (outfile)
print (outfile_bbox)

# s3 endpoint
#s3_base = "s3://lynker-spatial/"

##########################################################################
# Step 1: build subset
## caching the downloaded VPU files to "data" and writing all layers to "outfile"
div_file = subset_network(hl_uri = hl, cache_dir = "data", outfile = outfile)

# read and visualize the downloaded geopackage
div = read_sf(outfile, 'divides')
mapview(div)

# If you already have a basin geopackage
huc12_id_scan = '100302030705'
r_path = glue("/Users/ahmadjan/Core/simulations/owp/AGU_2023/HUC12-{huc12_id_scan}")
setwd(r_path)
outfile = glue('data/HUC12_{huc12_id_scan}.gpkg')
print (outfile)
div = read_sf(outfile, 'divides')
mapview(div)

##########################################################################
# Step 2: get DEM (maybe not needed in the future)
write_basin_dem <- function(path, gpkg_file, dem_path, dem_name) {
  elev_buf = read_sf(gpkg_file, 'divides') %>% st_union() %>% st_buffer(10000)
  elev_buf = get3DEP(elev_buf)
  plot(elev_buf$`30m CONUS DEM`)
  
  if (!dir.exists(dem_path)){
    dir.create(dem_path)
  }
  
  #name_dem = glue('{dem_name}.tif')
  writeRaster(elev_buf$`30m CONUS DEM`, filename=glue('{dem_path}/{dem_name}.tif'), overwrite=TRUE)  
}

dem_path = file.path(path,"DEM")
write_basin_dem(r_path, gpkg_file=outfile, dem_path, dem_name=gage_id)

##########################################################################
# Step 3: add cfe_noah attribtes to the geopackage
add_cfe_noah_attrs <- function(subset_path) {
  # Read all divide_id for the given hydrofabric subset
  divides <- sf::read_sf(subset_path, query = "SELECT divide_id FROM divides")$divide_id
  
  # Get all the CFE/NOAH-OWP attributes for the given catchments/divides
  vpu = unique(pull(read_sf(subset_path, "network")))
  #vpu = read_sf(outfile, "network") %>% pull(vpu) %>% unique()
  
  # full conus model param file is here: s3://lynker-spatial/v20/conus_model_attributes.parquet
  print (glue("VPU+: {vpu}"))
  #cfe <- arrow::open_dataset(glue("s3://nextgen-hydrofabric/pre-release/nextgen_{vpu}_cfe_noahowp.parquet")) |>
  #  dplyr::filter(divide_id %in% divides) |> 
  #  dplyr::collect()
  
  #cfe = arrow::read_parquet('s3://lynker-spatial/v20/model_attributes/nextgen_01.parquet') |>
  #  select(divide_id, twi_dist_4)
  
  cfe <- arrow::read_parquet(glue('s3://lynker-spatial/v20/model_attributes/nextgen_{vpu}.parquet')) |>
    dplyr::filter(divide_id %in% divides) |> 
    dplyr::collect()
  
  # Write the attributes to a new table in the hydrofabric subset GPKG
  sf::st_write(cfe, subset_path, layer = "model_attributes", append = FALSE)
}

# print layers before appending cfe attributes
layers_before_cfe_attr <- sf::st_layers(outfile)
print (layers_before_cfe_attr$name)

add_cfe_noah_attrs(outfile)

layers_after_cfe_attr <- sf::st_layers(outfile)
print (layers_after_cfe_attr$name)

##########################################################################
# Step 4: get topographic wetness index
add_twi <- function(subset_path) {
  
  twi <- dap("/vsis3/lynker-spatial/gridded-resources/twi.vrt", AOI = sf::read_sf(subset_path,'divides'))  
  # truncate negative values
  twi[twi < 0] <- 0
  twi[twi > 50] <- 50
  
  nclasses = 5
  twi_cat = execute_zonal(data = twi,
                          geom = div,
                          ID = "divide_id",
                          fun = zonal::equal_population_distribution,
                          groups = nclasses)
  print (twi_cat)
  twi_cat_breaks <<- execute_zonal(data = twi,
                                 geom = div,
                                 ID = "divide_id",
                                 fun = zonal::distribution,
                                 breaks = 30)
  #sf::st_write(twi_cat, outfile, layer = "twi", append = FALSE)
}

add_twi(outfile)


### TWI post-process/vis
#print (twi)
plot(twi)
plot(vect(outfile, "divides"), add = TRUE)

twi_vrt <- fromJSON(twi_cat_breaks["fun.twi"]$fun.twi[1])
print (twi_vrt)
#hist(r_list$frequency, breaks=5)

#freq = seq(1.0/nclasses,1, 1.0/nclasses)
barplot.default(rev(twi_vrt$v))#,  freq, names.arg = freq)



layers_after_twi_attr <- sf::st_layers(outfile)
print (layers_after_twi_attr$name)
