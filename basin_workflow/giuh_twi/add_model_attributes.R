
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
  
  model_attr <- arrow::read_parquet(glue('s3://lynker-spatial/v20/model_attributes/nextgen_{vpu}.parquet')) |>
    dplyr::filter(divide_id %in% divides) |> 
    dplyr::collect()
  
  # Write the attributes to a new table in the hydrofabric subset GPKG
  sf::st_write(model_attr, div_path, layer = "model_attributes", append = FALSE)
  
  return(model_attr)

}