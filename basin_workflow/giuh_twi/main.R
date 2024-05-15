# @author Ahmad Jan Khattak
# @email ahmad.jan@noaa.gov
# @date  December 22, 2023

# The script downloads geopackge(s) given USGS gauge id(s) (also can read gpkg from the disk)
# Computes TWI and GIUH and extracts model attributes using a parquet file on S3 endpoint
# Appends TWI, GIUH, and model attributes to .gpkg to be used for generating models' configuration files 

######################## REQUIRED INPUT #######################################
# Key steps are hihglight as:
# STEP #1: Setup (REQUIRED) (main.R)
# STEP #2: Options: Provide gage ID or list of gage IDs or set path to work with already download geopackages (main.R)
# Workflow substeps
#   STEP #3: Download geopackage (if needed) (inside driver.R)
#   STEP #4: Add model attributes to the geopackage (inside driver.R)
#   STEP #5: Compute TWI and width function (inside driver.R)
#   STEP #6: Compute GIUH (inside driver.R)
#   STEP #7: Compute Nash cascade parameters (N and K) for surface runoff (inside driver.R)
#   STEP #8: Append GIUH, TWI, width function, and Nash cascade parameters to model_attributes layer (inside driver.R)
###############################################################################


################################ SETUP #########################################
# STEP #1: INSTALL REQUIRED PACKAGES 
#          a) Install hydrofabric, zonal, etc. (if not already installed), and other
#          built-in libraries
#          b) Load custom .R files/functions such as giuh.R that computes GIUH
#          c) Specify DEM input file (Default points to S3 endpoint; see driver.R)
#          d) Set output directory path (output_dir)
################################################################################
# Point r_path to the directory of R scripts downloaded from the repository
r_path = "~/Core/SimulationsData/preprocessing/hydrofabric/smp_basin_workflow/basin_workflow/giuh_twi"
# (a)
source(paste0(r_path, "/install_load_libs.R"))
# (b)
source(glue("{r_path}/custom_functions.R"))
# (c)
#dem_infile = "/vsicurl/https://lynker-spatial.s3.amazonaws.com/gridded-resources/dem.vrt"

# (d) Point root_outpath to the directory where geopackage and other related files will be stored
#outpath_dir = "/Users/ahmadjan/Core/SimulationsData/preprocessing/CAMELS_2024/"
output_dir = "/Users/ahmadjan/Core/SimulationsData/preprocessing/test"
setwd(output_dir)
wbt_wd(getwd())

# create directory to stored catchment geopackage in case of errors or missing data
failed_dir = "failed_cats"
dir.create(failed_dir, recursive = TRUE, showWarnings = FALSE)


################################ OPTION #########################################
# STEP #2:Currently, the script works with a single gage ID, a list of gage ID, or already download
# geopackages (see Examples 1 and 2 below). 
# Once an option is selection, go to the corresponding Example and make sure to adjust the file/directories 
# to your local setting

###############################################################################
# Example 1: User-provided gage IDs; turn using_gage_IDs ON
option_using_gage_IDs <- TRUE
################################################################################
# Example 2: User-provided geopackages; turn using_gpkgs ON
option_using_gpkgs    <- FALSE



if (option_using_gage_IDs == TRUE) {
  ################################################################################
  # For this example either provide a gage ID or a file to read gage IDs from
  # Modify this part according your settings
  ################################################################################
  
  IDs_from_file <- TRUE
  gage_id <- NULL
  
  if (IDs_from_file) {
    
    gage_ids_infile = "/Users/ahmadjan/Core/SimulationsData/preprocessing/CAMELS_2024/CAMELS_v3_calib_BestModelsAll.csv"
    d = read.csv(gage_ids_infile,colClasses = c("character"))
    gage_ids <- d[["hru_id_CAMELS"]]    
    
  }
  else {
    gage_ids <- gage_id
  }
  
  stopifnot( length(gage_ids) > 0)
  
  cats_failed <- run_given_gage_IDs(gage_ids[c(1:2)], output_dir)      
  
  
} else if (option_using_gpkgs == TRUE) {
  ################################################################################
  # For this example set gpkg_i_dir to the directory containing geopackage(s)
  # Modify this part according your settings
  ################################################################################

  gpkgs_i_dir = glue("{output_dir}/camels_basins_all")   # input
  gage_files = list.files(gpkgs_i_dir, full.names = FALSE, pattern = "Gage_")
  
  cats_failed <- run_given_gpkg(gage_files, gpkgs_i_dir, output_dir)
  
}


print(cats_failed)

###############################################################################
# DONE
###############################################################################








#model_attr <- arrow::read_parquet(glue('s3://lynker-spatial/hydrofabric/v20.1/model_attributes/nextgen_01.parquet'))
#dplyr::as_tibble(model_attr)

# Check if the IDs exist in the HF conus network
#conus_net <- arrow::read_parquet(glue("{root_outpath}/model_attributes.parquet")) 
#|> 
#  dplyr::filter(divide_id %in% divides) |> 
#  dplyr::collect()


