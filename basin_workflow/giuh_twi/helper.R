# @author Ahmad Jan Khattak
# @email ahmad.jan@noaa.gov
# @date  February 05, 2024

# Get the DEM

dem_function <- function(div_infile,
                         dem_infile = "/vsicurl/https://lynker-spatial.s3.amazonaws.com/gridded-resources/dem.vrt",
                         dem_output_dir) {
  
  elev <- rast(dem_infile)
  
  # Get the catchment geopackage
  div <- read_sf(div_infile, 'divides')
  river <- read_sf(div_infile, "flowpaths")
  
  # Buffer because we want to guarantee we don not have boundary issues when processing the DEM
  div_bf <- st_buffer(div,dist=5000)
  
  dem <- crop(elev, project(vect(div_bf), crs(elev)), snap = "out")
  cm_to_m <- 0.01
  dem <- dem * cm_to_m
  writeRaster(dem, glue("{dem_output_dir}/dem.tif"), overwrite = TRUE)
  
  gdal_utils("warp",
             source = glue("{dem_output_dir}/dem.tif"),
             destination = glue("{dem_output_dir}/dem_proj.tif"),
             options = c("-of", "GTiff", "-t_srs", "EPSG:5070", "-r", "bilinear")
  )
  
  wbt_breach_depressions(dem = glue("{dem_output_dir}/dem_proj.tif"), output = glue("{dem_output_dir}/dem_corr.tif") )
  
}


#the condition [coverage_fraction > .1] excludes/drops all cell X that has fraction less than 10% in the divide Y
fun_crop_lower <- function(values, coverage_fraction) {
  data = (values * coverage_fraction)[coverage_fraction > 0.1]
  percentile_10 <- unname(quantile(data, probs = 0.15, na.rm = TRUE)) # unname function returns the quantile value only, and not the cut points
  data[data <= percentile_10] = percentile_10
}

fun_crop_upper <- function(values, coverage_fraction) {
  data = (values * coverage_fraction)[coverage_fraction > .1]
  percentile_90 <- unname(quantile(data, probs = 0.85, na.rm = TRUE))
  data[data >= percentile_90] = percentile_90
}


# Add model attribtes to the geopackage
add_model_attributes <- function(div_path) {
  # Read all divide_id for the given hydrofabric subset
  divides <- sf::read_sf(div_path, query = "SELECT divide_id FROM divides")$divide_id
  
  # Get all the CFE/NOAH-OWP attributes for the given catchments/divides
  vpu = unique(pull(read_sf(div_path, "network")))
  #vpu = read_sf(outfile, "network") %>% pull(vpu) %>% unique()
  
  # full conus model param file is here: s3://lynker-spatial/v20/conus_model_attributes.parquet
  print (glue("VPU+: {vpu}"))
  
  
  #cfe = arrow::read_parquet('s3://lynker-spatial/v20/model_attributes/nextgen_01.parquet') |>
  #  select(divide_id, twi_dist_4)
  hf_version = 'v20.1' 
  #arrow::read_parquet(glue('s3://lynker-spatial/v20/model_attributes/nextgen_{vpu}.parquet')) #older version v20
  
  model_attr <- arrow::read_parquet(glue('s3://lynker-spatial/hydrofabric/{hf_version}/model_attributes/nextgen_{vpu}.parquet')) |>
    dplyr::filter(divide_id %in% divides) |> 
    dplyr::collect()

  #if (nrow(model_attr) == 0) {
  #  model_attr <- arrow::read_parquet(glue('s3://lynker-spatial/{hf_version}/model_attributes.parquet')) |>
  #    dplyr::filter(divide_id %in% divides) |> 
  #    dplyr::collect()    
  #}

  #model_attr <- arrow::read_parquet(glue('s3://lynker-spatial/hydrofabric/{hf_version}/model_attributes.parquet')) |>
  #  dplyr::filter(divide_id %in% divides) |> 
  #  dplyr::collect() 
  
  #print ("Model attr B")
  #print (model_attr)
  cat ("m_attr: ", nrow(model_attr))
  stopifnot(nrow(model_attr) > 0)
  
  # Write the attributes to a new table in the hydrofabric subset GPKG
  sf::st_write(model_attr, div_path, layer = "model_attributes", append = FALSE)
  
  return(model_attr)
  
}