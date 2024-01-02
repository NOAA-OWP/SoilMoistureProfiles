

# Function computes topographic wetness index (TWI) and generates files 
# needed for computing width function (see below)
twi_function <- function(infile, directory, dem_path, nclasses = 5) {
  
  elev <- rast(dem_path)
  
  # Get the DEMph
  div <- read_sf(infile, 'divides')
  div_bf <- st_buffer(div,dist=5000)
  
  dem <- crop(elev, project(vect(div_bf), crs(elev)), snap = "out")
  cm_to_m <- 0.01
  dem <- dem * cm_to_m
  writeRaster(dem, glue("{directory}/dem.tif"), overwrite = TRUE)
  
  gdal_utils("warp",
             source = glue("{directory}/dem.tif"),
             destination = glue("{directory}/dem_proj.tif"),
             options = c("-of", "GTiff", "-t_srs", "EPSG:5070", "-r", "bilinear")
  )
  
  #dem_proj = raster(glue("{directory}/dem_proj.tif"))
  wbt_breach_depressions(dem = glue("{directory}/dem_proj.tif"), output = glue("{directory}/dem_corr.tif") )
  
  #################### TWI ########################
  
  # @param out_type Output type; one of 'cells' (default), 'catchment area', and 'specific contributing area'.
  wbt_d8_flow_accumulation(input = glue("{directory}/dem_corr.tif"), output = glue("{directory}/sca.tif")
                           , out_type = 'specific contributing area')
  
  
  wbt_slope(dem = glue("{directory}/dem_corr.tif"), output = glue("{directory}/slope.tif"))
  
  wbt_wetness_index(sca = glue("{directory}/sca.tif"), slope = glue("{directory}/slope.tif"), 
                    output = glue("{directory}/twi.tif"))
  
  
  twi = rast(glue("{directory}/twi.tif"))
  
  twi[twi < 0] <- 0
  twi[twi > 50] <- 50
  
  twi_cat <- execute_zonal(data = twi,
                            geom = div,
                            ID = "divide_id",
                            fun = zonal::equal_population_distribution,
                            groups = nclasses)
  #print (twi_cat)
  
  #twi_compute <<- execute_zonal(data = twi,
  #                              geom = div,
  #                              ID = "divide_id",
  #                              fun = zonal::distribution,
  #                              breaks = 30)
  

  return(twi_cat)
}



width_function <- function(infile, directory) {
  
  div <- read_sf(infile, 'divides')
  
  wbt_d8_pointer(dem = glue("{directory}/dem_corr.tif"), 
                 output = glue("{directory}/dem_d8.tif"))
  
  wbt_downslope_flowpath_length(d8_pntr = glue("{directory}/dem_d8.tif"), 
                                output = glue("{directory}/down_fp_length.tif"), watersheds=NULL)
  # note: watersheds=div never tested but maybe useful in some cases; defualt is NULL
  
  #dem_sca <- raster(glue("{directory}/sca.tif"))
  flowpath_length <- rast(glue("{directory}/down_fp_length.tif"))
  
  fp_min_ftn = execute_zonal(data = flowpath_length,
                             geom = div,
                             ID = "divide_id",
                             fun = fun_crop_lower)
  
  
  fp_max_ftn = execute_zonal(data = flowpath_length,
                             geom = div,
                             ID = "divide_id",
                             fun = fun_crop_upper)  
  
  # create a grid based on min values of sub-catchments, grid resolution consistent with downslope_flowpth_length
  rasterized_fp_min <- rasterize(fp_min_ftn, flowpath_length, field=fp_min_ftn$fun.down_fp_length)
  rasterized_fp_max <- rasterize(fp_max_ftn, flowpath_length, field=fp_max_ftn$fun.down_fp_length)
  
  # assign fp_length raster to rast_fp_temp 
  rast_fp_temp <- flowpath_length
  
  # assign NA to values (per catchment) less than the min and max
  rast_fp_temp[rast_fp_temp <= rasterized_fp_min] <- NA
  rast_fp_temp[rast_fp_temp >= rasterized_fp_max] <- NA
  
  # subtract adjusted-raster (created 15% and 85% quantiles) from minimun-catchment-value to get localized distance
  # of a point to the catchment outlet
  downslope_fp_distance_cat_outlet <- rast_fp_temp - rasterized_fp_min
  
  # channel cumulative distribution of area with distance
  width_dist <- execute_zonal(data = downslope_fp_distance_cat_outlet,
                              geom = div,
                              ID = "divide_id",
                              fun = zonal::distribution,
                              breaks = c(0.0,0.1,500,1000,1500))
  
  return(width_dist)
}

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


### Pre-computed TWI values provided in the hydrofabric
twi_width_pre_computed <- function(div_path) {
  
  twi <- dap("/vsis3/lynker-spatial/gridded-resources/twi.vrt", AOI = sf::read_sf(div_path,'divides'))  
  # truncate negative values
  twi[twi < 0] <- 0
  twi[twi > 50] <- 50
  
  twi_cat = execute_zonal(data = twi,
                          geom = div,
                          ID = "divide_id",
                          fun = zonal::equal_population_distribution,
                          groups = nclasses)
  return(twi_cat)
}

