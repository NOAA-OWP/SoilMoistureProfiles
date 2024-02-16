# @author Ahmad Jan Khattak
# @email ahmad.jan@noaa.gov
# @date  February 05, 2024

# ########################### TWI ########################
# Function computes topographic wetness index (TWI) and generates files 
# needed for computing width function (see below)
twi_function <- function(infile, directory, distribution = 'quantiles', nclasses = 5) {
  
  div <- read_sf(infile, 'divides')
  
  # @param out_type Output type; one of 'cells' (default), 'catchment area', and 'specific contributing area'.
  wbt_d8_flow_accumulation(input = glue("{directory}/dem_corr.tif"), output = glue("{directory}/sca.tif")
                           , out_type = 'specific contributing area')
  
  
  wbt_slope(dem = glue("{directory}/dem_corr.tif"), output = glue("{directory}/slope.tif"))
  
  wbt_wetness_index(sca = glue("{directory}/sca.tif"), slope = glue("{directory}/slope.tif"), 
                    output = glue("{directory}/twi.tif"))
  
  twi = rast(glue("{directory}/twi.tif"))
  print (twi)
  twi[twi < 0] <- 0
  twi[twi > 50] <- 50
  
  if (distribution == 'quantiles') {
    twi_cat <- execute_zonal(data = twi,
                             geom = div,
                             ID = "divide_id",
                             fun = equal_population_distribution,
                             groups = nclasses)
    
  }
  else if (distribution == 'simple') {
     twi_cat <- execute_zonal(data = twi,
                              geom = div,
                              ID = "divide_id",
                              fun = zonal::distribution,
                              breaks = nclasses)
  }

  return(twi_cat)
}



width_function <- function(infile, directory) {
  
  div <- read_sf(infile, 'divides')
  
  wbt_d8_pointer(dem = glue("{directory}/dem_corr.tif"), 
                 output = glue("{directory}/dem_d8.tif"))
  
  wbt_downslope_flowpath_length(d8_pntr = glue("{directory}/dem_d8.tif"), 
                                output = glue("{directory}/downslope_fp_length.tif"), watersheds=NULL)
  # note: watersheds=div never tested but maybe useful in some cases; defualt is NULL
  
  flowpath_length <- rast(glue("{directory}/downslope_fp_length.tif"))
  
  fp_min_ftn = execute_zonal(data = flowpath_length,
                             geom = div,
                             ID = "divide_id",
                             fun = fun_crop_lower)
  
  
  fp_max_ftn = execute_zonal(data = flowpath_length,
                             geom = div,
                             ID = "divide_id",
                             fun = fun_crop_upper)  
  
  # create a grid based on min values of sub-catchments, grid resolution consistent with downslope_flowpth_length
  rasterized_fp_min <- rasterize(fp_min_ftn, flowpath_length, field=fp_min_ftn$fun.downslope_fp_length)
  rasterized_fp_max <- rasterize(fp_max_ftn, flowpath_length, field=fp_max_ftn$fun.downslope_fp_length)
  
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

