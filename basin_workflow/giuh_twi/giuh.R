# @author Ahmad Jan Khattak
# @email ahmad.jan@noaa.gov
# @date  February 05, 2024

# NOTE: # TauDem uses Dinf:  D-infinity flow direction is defined as steepest downward slope 
#          on planar triangular facets on a block centered grid. Here we use D8 scheme


# Function computes Geomorphological Instantaneous Unit Hydrograph (GIUH)
giuh_function <- function(div_infile, dem_output_dir, vel_channel = 1, vel_overland = .5, vel_gully = .2, gully_threshold = 3) {
  
  div <- read_sf(div_infile, 'divides')
  river <- read_sf(div_infile, "flowlines")
  
  # @param out_type Output type; one of 'cells' (default), 'catchment area', and 'specific contributing area'.
  wbt_d8_flow_accumulation(input = glue("{dem_output_dir}/dem_corr.tif"), output = glue("{dem_output_dir}/giuh_sca.tif"),
                           out_type = 'specific contributing area')
  
  sca <- rast(glue("{dem_output_dir}/giuh_sca.tif"))
  rasterized_river <- rasterizeGeom(vect(river), sca, fun="length")
  
  writeRaster(rasterized_river, glue("{dem_output_dir}/giuh_river.tif"), overwrite = TRUE)
  
  #x <- ifel(sca <= gully_threshold, vel_gully, vel_overland) #original script
  x <- ifel(sca > gully_threshold, vel_gully, vel_overland)
  x <- ifel(rasterized_river > 0, vel_channel, x)
  
  sec_to_min = 60  # 60 minutes (1 hour) is the time discretization for giuh
  # This x matrix represents the time (in minutes) it takes to travel 1 meter. 
  # Speed = distance[meters]/Time[minutes]. X is 1/Speed 
  x <- 1.0/(x*sec_to_min) 
  names(x)  = "travel_time"
  
  writeRaster(x, glue("{dem_output_dir}/giuh_travel_time.tif") ,overwrite=TRUE)  
  
  # This one calculates the path to the basin outlet
  wbt_d8_pointer(dem = glue("{dem_output_dir}/dem_corr.tif"), output = glue("{dem_output_dir}/dem_d8.tif"))
  
  # Using S = V * T => T = S/V; divide distance (flowpath_length) by weights (1/V)
  wbt_downslope_flowpath_length(d8_pntr = glue("{dem_output_dir}/dem_d8.tif"),
                                output  = glue("{dem_output_dir}/giuh_minute.tif"),
                                weights = glue("{dem_output_dir}/giuh_travel_time.tif"))
  
  # from basin outlet to catchment outlet workflow
  giuh_minute <- rast(glue("{dem_output_dir}/giuh_minute.tif"))
  
  time_min_ftn = execute_zonal(data = giuh_minute,
                               geom = div,
                               ID = "divide_id",
                               fun = fun_crop_lower)
  
  
  time_max_ftn = execute_zonal(data = giuh_minute,
                               geom = div,
                               ID = "divide_id",
                               fun = fun_crop_upper)  
  
  # create a grid based on min values of sub-catchments, grid resolution consistent with downslope_flowpth_length
  rasterized_time_min <- rasterize(time_min_ftn, giuh_minute, field=time_min_ftn$fun.giuh_minute)
  rasterized_time_max <- rasterize(time_max_ftn, giuh_minute, field=time_max_ftn$fun.giuh_minute)
  
  # assign giuh_minute raster to rast_giuh_temp 
  rast_giuh_minute_temp <- giuh_minute
  
  # assign NA to values (per catchment) less than the min and max
  rast_giuh_minute_temp[rast_giuh_minute_temp <= rasterized_time_min] <- NA
  rast_giuh_minute_temp[rast_giuh_minute_temp >= rasterized_time_max] <- NA
  
  # subtract adjusted-raster (created 15% and 85% quantiles) from minimun-catchment-value to get localized distance
  # of a point to the catchment outlet
  downslope_giuh_cat_outlet <- rast_giuh_minute_temp - rasterized_time_min
  
  #print (downslope_giuh_cat_outlet)
  
  writeRaster(downslope_giuh_cat_outlet, glue("{dem_output_dir}/downslope_giuh_cat_outlet.tif"),
              overwrite=TRUE)  
  
  # channel cumulative distribution of area with distance
  #giuh_dist <- execute_zonal(data = downslope_giuh_cat_outlet,
  #                           geom = div,
  #                           ID = "divide_id",
  #                           fun = zonal::distribution,
  #                           breaks = seq(0.0, 600, by=60),
  #                           constrain = TRUE)
  
  giuh_dist <- execute_zonal(data = downslope_giuh_cat_outlet,
                             geom = div,
                             ID = "divide_id",
                             fun = dist2,
                             breaks = seq(0.0, 600, by=60),
                             constrain = TRUE)
  
  return(giuh_dist)
}

dist2 <- function(value, coverage_fraction, breaks = 10, constrain = FALSE) {
  if (length(value) <= 0 | all(is.nan(value))) {
    return("[]")
  }
  
  x1 = value * coverage_fraction
  x1 = x1[!is.na(x1)]
  
  if (constrain | length(breaks) > 1) {
    breaks_tmp = c(breaks[1],breaks[2])
    #message(paste0("BREAKS PRIOR: ", paste(breaks, collapse = ", ")))
    breaks = breaks[breaks <= max(x1, na.rm = TRUE)]
    if (length(breaks) == 1) {
      breaks = breaks_tmp
    }
    #message(paste0("BREAKS AFTER: ", paste(breaks, collapse = ", ")))
  }
  
  
  tmp = as.data.frame(table(cut(x1, breaks = breaks)))
  
  tmp$v = as.numeric(gsub("]", "", sub('.*,\\s*', '', tmp$Var1)))
  
  tmp$frequency = tmp$Freq / sum(tmp$Freq)
  
  as.character(toJSON(tmp[,c("v", "frequency")]))
}




