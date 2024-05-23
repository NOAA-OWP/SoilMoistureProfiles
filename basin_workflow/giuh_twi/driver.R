######################### DOWNLOAD GEOPACKAGE ##################################
# STEP #3: provide USGS gauge id or your own geopackage (single or multiple)
################################################################################

run_driver <- function(gage_id = NULL, is_gpkg_provided = FALSE, 
                       dem_infile = "/vsicurl/https://lynker-spatial.s3.amazonaws.com/gridded-resources/dem.vrt", 
                       dem_output_dir,
                       loc_gpkg_file = "") {
  
  outfile <- " "
  if(!is_gpkg_provided) {
    #gage_id = '01033000' # A hydrolocation URI
    #hl = glue('Gages-{gage_id}')
    fid = glue('USGS-{gage_id}')
    outfile <- glue('data/gage_{gage_id}.gpkg')
    #cat ("run_main:", fid, outfile, "\n")
    ## caching the downloaded VPU files to "data" and writing all layers to "outfile"
    #subset_network(hl_uri = hl, cache_dir = "data", outfile = outfile)
    
    subset_network(nldi_feature = list(featureSource="nwissite", featureID=fid),
                   cache_dir = "data", outfile = outfile,
                   base_s3 = "s3://lynker-spatial/hydrofabric/v20.1")
  } else {
    outfile <- loc_gpkg_file
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
  # Note: The default distribution = 'quantiles'
  ###############################################################################

  dem_function(div_infile = outfile, dem_infile, dem_output_dir)
  
  twi <- twi_function(div_infile = outfile, dem_output_dir = dem_output_dir, distribution = 'simple', nclasses = 30)
  
  width_dist <- width_function(div_infile = outfile, dem_output_dir = dem_output_dir)
  
  twi_dat_values = data.frame(ID = twi$divide_id, twi = twi$fun.twi, width_dist = width_dist$fun.downslope_fp_length)
  
  # write TWI and width function layers to the geopackage
  names(twi_dat_values)
  colnames(twi_dat_values) <- c('divide_id', 'twi', 'width_dist')
  names(twi_dat_values)
  #sf::st_write(dat_values, outfile, layer = "twi", append = FALSE)
  
  ### NOTES: Pre-computed TWI
  # Note 1: model attributes layer ships with pre-computed TWI distribution with four equal quantiles
  #m_attr$twi_dist_4
  
  # Note 2: The user can also compute their own distribution from the pre-computed TWI using the dataset
  # available at s3://lynker-spatial/gridded-resources/twi.vrt
  
  twi_pre_computed <- twi_pre_computed_function(div_infile = outfile, distribution = 'simple', nclasses = 30)
  
  ############################### GENERATE GIUH ##################################
  # STEP #6: Generate GIUH and write to the geopackage
  ################################################################################
  # There are many "model" options to specify the velocity.
  # Here we are using a simple approach: constant velocity as a function of upstream drainage area.
  vel_channel     <- 1.0  # meter/second
  vel_overland    <- 0.1  # Fred: 0.1
  vel_gully       <- 0.2 # meter per second
  gully_threshold <- 30.0 # m (longest , closer to 10-30 m, Refs) 
  
  giuh_compute <- giuh_function(div_infile = outfile, dem_output_dir = dem_output_dir, 
                                vel_channel, vel_overland, vel_gully, gully_threshold)
  
  #giuh_compute[2,] %>% t()
  
  # write GIUH layer to the geopackage
  giuh_dat_values = data.frame(ID = giuh_compute$divide_id, giuh = giuh_compute$fun.giuh_minute)
  names(giuh_dat_values)
  colnames(giuh_dat_values) <- c('divide_id', 'giuh')
  names(giuh_dat_values)
  
  #giuh_dat_values$giuh[1]
  
  
  ############################### GENERATE GIUH ##################################
  # STEP #7: Generate Nash cascade parameters for surface runoff
  ################################################################################
  nash_params_surface <- get_nash_params(giuh_dat_values, calib_n_k = FALSE)
  
  
  ################################################################################
  # STEP #8: Append GIUH, TWI, width function, and Nash cascade N and K parameters
  # to model attributes layers
  ################################################################################
  m_attr$giuh <- giuh_dat_values$giuh # append GIUH column to the model attributes layer
  m_attr$twi <- twi_dat_values$twi   # append TWI column to the model attributes layer
  m_attr$width_dist <- twi_dat_values$width_dist # append width distribution column to the model attributes layer
  
  m_attr$N_nash_surface <- nash_params_surface$N_nash
  m_attr$K_nash_surface <- nash_params_surface$K_nash
  
  sf::st_write(m_attr, outfile,layer = "model_attributes", append = FALSE)
  
}


###############################################################################
# main script that loops over all the gage IDs and computes giuh/twi etc.
run_given_gage_IDs <- function(gage_ids, output_dir) {
  # vector contains ID of basins that failed for some reason
  cats_failed <- numeric(0)
  
  for (id in gage_ids) {
    cat_dir = glue("{output_dir}/{id}")
    dir.create(cat_dir, recursive = TRUE, showWarnings = FALSE)
    setwd(cat_dir)
    wbt_wd(getwd())
    # DEM and related files (such as projected/corrected DEMs, and specific contributing area rasters are stored here)
    dem_o_dir = "dem"
    dir.create(dem_o_dir, recursive = TRUE, showWarnings = FALSE)
    
    failed <- TRUE
    tryCatch({
      cat ("Running ", id, "\n")
      run_driver(gage_id = id, dem_output_dir = dem_o_dir)
      
      failed <- FALSE
    }, error = function(e) {
      failed <- TRUE
    })
    
    if (failed) {
      cat ("Cat failed:", id, "\n")
      cats_failed <- append(cats_failed, id)
      file.rename(cat_dir, glue("{output_dir}/{failed_dir}/{id}"))
      next
    }
    else {
      cat ("Cat passed:", id, "\n")
      #cats_passed <- append(cats_passed, id)
    }
    
  }
  
  setwd(output_dir)
  
  return(cats_failed)
  
}

###############################################################################
# main script that loops over all the geopackages and computes giuh/twi etc.
run_given_gpkg <- function(gages_files, gpkgs_i_dir, output_dir) {

  # vector contains ID of basins that failed for some reason
  cats_failed <- numeric(0)
  
  print(gage_files)
  for (f in gage_files) {
    
    id = as.integer(unlist(strsplit(f, split = '[_&.]+'))[2])
    #cat (nchar(id), as.integer(nchar(id)/2) , as.integer(nchar(id)/2) %% 2, "\n")
    # check if gage ID is missing a leading zero, does not happens most of the times, but good to check
    if (as.integer(nchar(id)/2) %% 2 == 1) {
      id <- paste(0,id, sep = "")
    }
    
    #if (nchar(id) < 8) {
    #  id <- paste(0,id, sep = "")
    #}
    
    print (id)


    cat_o_dir = glue("{output_dir}/{id}")
    print (cat_o_dir)
    
    dir.create(cat_o_dir, recursive = TRUE, showWarnings = FALSE)
    setwd(cat_o_dir)
    wbt_wd(getwd())
    
    # DEM and related files (such as projected/corrected DEMs, and specific contributing 
    # area rasters are stored here)
    dem_o_dir = "dem"
    dir.create(dem_o_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create("data", recursive = TRUE, showWarnings = FALSE)
    file.copy(glue("{gpkgs_i_dir}/{f}"), "data")
    file.rename(glue("data/{f}"), glue("data/Gage_{id}.gpkg"))
    
    failed <- TRUE
    
    tryCatch({
      cat ("Running ", id, "\n")
      print (glue("{cat_o_dir}/data/Gage_{id}.gpkg"))
      run_driver(is_gpkg_provided = TRUE, loc_gpkg_file = glue("{cat_o_dir}/data/Gage_{id}.gpkg"), 
                 dem_output_dir = dem_o_dir)
      
      failed <- FALSE
      
    }, error = function(e) {
      failed <- TRUE
    })
    
    if (failed) {
      cat ("Cat failed:", id, "\n")
      cats_failed <- append(cats_failed, id)
      file.rename(cat_o_dir, glue("{output_dir}/{failed_dir}/{id}"))
      next
    }
    else {
      cat ("Cat passed:", id, "\n")
      #cats_passed <- append(cats_passed, id)
    }
    
  } 
  
  setwd(output_dir)
  
  return(cats_failed)
}