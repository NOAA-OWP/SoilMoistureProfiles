############################################################################################
# Author  : Ahmad Jan
# Contact : ahmad.jan@noaa.gov
# Date    : September 28, 2023
############################################################################################

"""
The script generates config files for all nextgen models using a hydrofabric subset geopackage
 - inputs:  see main function for inputs (taken as arguments)
 - outputs: writes configuration files for all catchments within the basin
"""

import os, sys
import argparse
import re
import copy
import glob
import json
import subprocess
import pandas as pd
import geopandas as gpd
import numpy as np
import fiona

#############################################################################
# module reads NWM soil type file and returns a table
# this is used for two purposes: 1) Xinanjiang parameters, 2) soil quartz content used in soil freeze thaw model
# @param infile : input file contain NWM 3.0 soil data/properties
# - returns     : dataframe
#############################################################################
def get_soil_class_NWM(infile):
    header = ['index','BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI','SATDK','SATDW','WLTSMC', \
              'QTZ', 'BVIC', 'AXAJ', 'BXAJ', 'XXAJ', 'BDVIC', 'BBVIC', 'GDVIC','ISLTYP']
    
    df = pd.read_table(infile, delimiter=',', index_col=0, skiprows=3, nrows=19, names=header)

    return df


#############################################################################
# module reads hydrofabric geopackage file and retuns a dict containing parameters needed for our models
# this is intended to be modified if more models are added or more soil parameters need to be extracted
# @param infile : input file pointing to hydrofabric basin geopkacge
# - returns     : geodataframe 
#############################################################################
def read_gpkg_file(infile, coupled_models):

    gdf_soil = gpd.read_file(infile, layer='model_attributes')
    gdf_soil.set_index("divide_id", inplace=True)

    layers = fiona.listlayers(infile)
    print ("Geopackage layers: ", layers)
    print ("\n")
    
    # extract TWI for topmodel
    if ("nom_topmodel" in coupled_models):
        gdf_twi = gpd.read_file(infile, layer='twi')
        gdf_twi.set_index("divide_id", inplace=True)

    if ("cfe" in coupled_models or "lasam" in coupled_models):
        gdf_giuh = gpd.read_file(infile, layer='giuh')
        gdf_giuh.set_index("divide_id", inplace=True)
        
    gdf_soil['bexp_soil_layers_stag=1'].fillna(16,inplace=True)
    gdf_soil['dksat_soil_layers_stag=1'].fillna(0.00000338,inplace=True)
    gdf_soil['psisat_soil_layers_stag=1'].fillna(0.355,inplace=True)  
    gdf_soil['smcmax_soil_layers_stag=1'].fillna(0.439,inplace=True)
    gdf_soil['smcwlt_soil_layers_stag=1'].fillna(0.066,inplace=True)
    gdf_soil['gw_Zmax'].fillna(0.01,inplace=True)
    gdf_soil['gw_Coeff'].fillna(1.8e-05,inplace=True)
    gdf_soil['gw_Expon'].fillna(6.0,inplace=True)
    gdf_soil['slope'].fillna(1.0, inplace=True)
    gdf_soil['ISLTYP'].fillna(1,inplace=True)
    gdf_soil['IVGTYP'].fillna(1,inplace=True)
    
    gdf_soil['gw_Zmax'] = gdf_soil['gw_Zmax']/1000. 
    gdf_soil['gw_Coeff'] = gdf_soil['gw_Coeff']*3600*pow(10,-6)
    
    
    if('refkdt' in gdf_soil):
        gdf_soil['refkdt'].fillna(3.0,inplace=True)
    else:
        gdf_soil['refkdt']  = 3.0

    
    # copy parameters needed
    gdf = gpd.GeoDataFrame(pd.DataFrame(), geometry= gdf_soil['geometry'], index=gdf_soil.index)
    gdf['soil_params.b']       = gdf_soil['bexp_soil_layers_stag=1'].copy()
    gdf['soil_params.satdk']   = gdf_soil['dksat_soil_layers_stag=1'].copy()
    gdf['soil_params.satpsi']  = gdf_soil['psisat_soil_layers_stag=1'].copy()
    gdf['soil_params.slop']    = gdf_soil['slope'].copy()
    gdf['soil_params.smcmax']  = gdf_soil['smcmax_soil_layers_stag=1'].copy()
    gdf['soil_params.wltsmc']  = gdf_soil['smcwlt_soil_layers_stag=1'].copy()
    gdf['refkdt']              = gdf_soil['refkdt'].copy()
    gdf['max_gw_storage']      = gdf_soil['gw_Zmax'].copy()
    gdf['Cgw']                 = gdf_soil['gw_Coeff'].copy()
    gdf['expon']               = gdf_soil['gw_Expon'].copy()
    gdf['ISLTYP']              = gdf_soil['ISLTYP'].copy()
    gdf['IVGTYP']              = gdf_soil['IVGTYP'].copy()

    # ensure parameter `b` is non-zero
    mask = gdf['soil_params.b'].gt(0.0) # greater than or equal to
    min_value = gdf['soil_params.b'][mask].min() # get the min value > 0.0

    mask = gdf['soil_params.b'].le(0.0) # find all values <= 0.0
    
    #df['soil_params.b'][mask] = min_value
    gdf.loc[mask, 'soil_params.b'] = min_value

    # TWI for topmodel
    if ("nom_topmodel" in coupled_models):
        gdf['twi'] = gdf_twi['twi']
        gdf['width_dist'] = gdf_twi['width_dist']
        
    if ("cfe" in coupled_models or "lasam" in coupled_models):
        gdf['giuh'] = gdf_giuh['giuh']
        
    # get catchment ids
    df_cats = gpd.read_file(infile, layer='divides')
    catids = [int(re.findall('[0-9]+',s)[0]) for s in df_cats['divide_id']]
    
    return gdf, catids

#############################################################################
# write uniform forcing data files, takes a base file and replicate it over all catchments (will be removed later)
# this is only for testing purposes
# @param catids : array/list of integers contain catchment ids
# @param infile : input file of the  base forcing data     
#############################################################################
def write_forcing_files(catids, infile):
    
    inpath= os.path.dirname(infile)
    
    for catID in catids:
        cat_name = 'cat-'+str(catID)
        fname = cat_name+'.csv'
        str_sub="cp " + infile + " "+os.path.join(inpath,fname)
        out=subprocess.call(str_sub,shell=True)

#############################################################################
# The function generates configuration files for NOAH-OWP-Modular

# @param catids      : array/list of integers contain catchment ids
# @param nom_dir     : output directory (config files are written to this directory)
# @param forcing_dir : forcing data directory containing data for each catchment
# @param gpkg_file   : basin geopackage file
# @param simulation_time : dictionary contain start/end time of the simulation

#############################################################################
def write_nom_input_files(catids, nom_dir, forcing_dir, gpkg_file, simulation_time):
    
    df_soil = gpd.read_file(gpkg_file, layer='model_attributes')
    df_cats = gpd.read_file(gpkg_file, layer='divides')
    df_cats = df_cats.to_crs("EPSG:4326") # change CRS to 4326
    
    df_soil.set_index("divide_id", inplace=True)
    df_cats.set_index("divide_id", inplace=True)

    # get soil type and fill with 1 if nan
    df_soil['ISLTYP'].fillna(1,inplace=True)

    print ("NOM simulation time: ", simulation_time)
    start_time = pd.Timestamp(simulation_time['start_time']).strftime("%Y%m%d%H%M")
    end_time   = pd.Timestamp(simulation_time['end_time']).strftime("%Y%m%d%H%M")
    
    for catID in catids:
        cat_name = 'cat-'+str(catID)
        
        centroid_x = str(df_cats['geometry'][cat_name].centroid.x)
        centroid_y = str(df_cats['geometry'][cat_name].centroid.y)

        soil_type = str(df_soil.loc[cat_name]['ISLTYP'])
        veg_type = str(df_soil.loc[cat_name]['IVGTYP'])
        
        timing = ["&timing                                   ! and input/output paths",
                  "  dt                 = 3600.0             ! timestep [seconds]",
                  "  startdate          = \"%s\"             ! UTC time start of simulation (YYYYMMDDhhmm)"%start_time,
                  "  enddate            = \"%s\"             ! UTC time end of simulation (YYYYMMDDhhmm)"%end_time,
                  "  forcing_filename   = \"%s.csv\"         ! file containing forcing data"%(os.path.join(forcing_dir,cat_name)),
                  "  output_filename    = \"output-%s.csv\""%cat_name,
                  "/\n"
                  ]
        
        params = ["&parameters",
                  "  parameter_dir      = \"%s\"  ! location of input parameter files"%(os.path.join(nom_dir,"parameters")),
                  "  general_table      = \"GENPARM.TBL\"                ! general param tables and misc params",
                  "  soil_table         = \"SOILPARM.TBL\"               ! soil param table",
                  "  noahowp_table      = \"MPTABLE.TBL\"                ! model param tables (includes veg)",
                  "  soil_class_name    = \"STAS\"                       ! soil class data source - STAS or STAS-RUC",
                  "  veg_class_name     = \"MODIFIED_IGBP_MODIS_NOAH\"   ! vegetation class data source - MODIFIED_IGBP_MODIS_NOAH or USGS",
                  "/\n"
                  ]
        
        location = ["&location                                         ! for point runs",
                    "  lat              = %s                           ! latitude [degrees]  (-90 to 90)"%centroid_y,
                    "  lon              = %s                           ! longitude [degrees] (-180 to 180)"%centroid_x,
                    "  terrain_slope    = 0.0                          ! terrain slope [degrees]",
                    "  azimuth          = 0.0                          ! terrain azimuth or aspect [degrees clockwise from north]",
                    "/ \n"
                    ]
        
        forcing = ["&forcing",
                   "  ZREF               = 10.0                        ! measurement height for wind speed (m)",
                   "  rain_snow_thresh   = 1.0                         ! rain-snow temperature threshold (degrees Celcius)",
                   "/ \n"
                   ]

        model_opt = ["&model_options                                   ! see OptionsType.f90 for details",
                     "  precip_phase_option               = 1",
                     "  snow_albedo_option                = 1",
                     "  dynamic_veg_option                = 4",
                     "  runoff_option                     = 3",
                     "  drainage_option                   = 8",
                     "  frozen_soil_option                = 1",
                     "  dynamic_vic_option                = 1",
                     "  radiative_transfer_option         = 3",
                     "  sfc_drag_coeff_option             = 1",
                     "  canopy_stom_resist_option         = 1",
                     "  crop_model_option                 = 0",
                     "  snowsoil_temp_time_option         = 3",
                     "  soil_temp_boundary_option         = 2",
                     "  supercooled_water_option          = 1",
                     "  stomatal_resistance_option        = 1",
                     "  evap_srfc_resistance_option       = 4",
                     "  subsurface_option                 = 2",
                     "/\n",
                     ]

        struct = ["&structure",
                  "  isltyp           = %s               ! soil texture class"%soil_type,
                  "  nsoil            = 4               ! number of soil levels",
                  "  nsnow            = 3               ! number of snow levels",
                  "  nveg             = 27              ! number of vegetation types",
                  "  vegtyp           = %s               ! vegetation type"%veg_type,
                  "  croptype         = 0               ! crop type (0 = no crops; this option is currently inactive)",
                  "  sfctyp           = 1               ! land surface type, 1:soil, 2:lake",
                  "  soilcolor       = 4               ! soil color code",
                  "/\n"
                  ]
        
        init_val = ["&initial_values",
                    "  dzsnso    =  0.0,  0.0,  0.0,  0.1,  0.3,  0.6,  1.0     ! level thickness [m]",
                    "  sice      =  0.0,  0.0,  0.0,  0.0                       ! initial soil ice profile [m3/m3]",
                    "  sh2o      =  0.3,  0.3,  0.3,  0.3                       ! initial soil liquid profile [m3/m3]",
                    "  zwt       =  -2.0                                        ! initial water table depth below surface [m]",
                    "/\n",
                    ]

        # combine all sub-blocks
        nom_params = timing + params + location + forcing + model_opt + struct + init_val

        fname_nom = cat_name+'_config_nom.input'
        nom_file = os.path.join(nom_dir, fname_nom)
        with open(nom_file, "w") as f:
            f.writelines('\n'.join(nom_params))
            

#############################################################################
# The function generates configuration file for CFE
# @param catids         : array/list of integers contain catchment ids
# @param runoff_schame  : surface runoff schemes - Options = Schaake or Xinanjiang
# @soil_class_NWM       : a dict containing NWM soil characteristics, used here for extracting
#                         Xinanjiang properties for a given soil type
# @gdf_soil             : geodataframe contains soil properties extracted from the hydrofabric
# @param cfe_dir        : output directory (config files are written to this directory)
# @param giuh_dir       : GIUH data directory (pre-computed GIUH distributions for each catchment)
# @param gpkg_file      : basin geopackage file
# @param coupled_models : option needed to modify CFE config files based on the coupling type
#############################################################################
def write_cfe_input_files(catids, runoff_scheme, soil_class_NWM, giuh_dir,
                          gdf_soil, cfe_dir, coupled_models):

    if (not runoff_scheme in ["Schaake", "Xinanjiang"]):
        sys.exit("Runoff scheme should be: Schaake or Xinanjiang")
    
    ice_content_threshold  = 0.3 # used when coupled with Soil freeze thaw model

    urban_decimal_fraction = 0.0 # used when runoff scheme is Xinanjiang
    
    # loop over all catchments and write config files
    for catID in catids:
        cat_name = 'cat-'+str(catID) 
        fname = cat_name+'*.txt'

        #giuh_file = glob.glob(os.path.join(giuh_dir, fname))[0]
        #df_giuh = pd.read_table(giuh_file,  delimiter='=', names=["Params","Values"], index_col=0)

        # cfe params set
        cfe_params = ['forcing_file=BMI',
                      'surface_partitioning_scheme=Schaake',
                      'soil_params.depth=2.0[m]',
                      'soil_params.b=' + str(gdf_soil['soil_params.b'][cat_name])+'[]',
                      'soil_params.satdk=' + str(gdf_soil['soil_params.satdk'][cat_name])+'[m s-1]', 
                      'soil_params.satpsi=' + str(gdf_soil['soil_params.satpsi'][cat_name])+'[m]',
                      'soil_params.slop='+str(gdf_soil['soil_params.slop'][cat_name])+"[m/m]",
                      'soil_params.smcmax='+str(gdf_soil['soil_params.smcmax'][cat_name])+"[m/m]",
                      'soil_params.wltsmc='+str(gdf_soil['soil_params.wltsmc'][cat_name])+"[m/m]",
                      'soil_params.expon=1.0[]',
                      'soil_params.expon_secondary=1.0[]',
                      'refkdt='+str(gdf_soil['refkdt'][cat_name]),
                      'max_gw_storage='+str(gdf_soil['max_gw_storage'][cat_name])+'[m]',
                      'Cgw='+str(gdf_soil['Cgw'][cat_name])+'[m h-1]',
                      'expon='+str(gdf_soil['expon'][cat_name])+'[]',
                      'gw_storage=0.05[m/m]',
                      'alpha_fc=0.33',
                      'soil_storage=0.1[m/m]',
                      'K_nash=0.03[]',
                      'K_lf=0.01[]',
                      'nash_storage=0.0,0.0',
                      'num_timesteps=1',
                      'verbosity=0',
                      'DEBUG=0'
                   ]

        
        if (gdf_soil['soil_params.b'][cat_name] == 1.0):
            cfe_list[3] = 1.1
        
        #cfe_lst[24] += df_giuh.loc['giuh_ordinates'].iloc[0]
        #cfe_params.append('giuh_ordinates='+df_giuh.loc['giuh_ordinates'].iloc[0])
        #cfe_params.append('giuh_ordinates=0.3,0.25,0.2,0.15,0.1')

        # add giuh ordinates
        giuh_cat = json.loads(gdf_soil['giuh'][cat_name])
        giuh_cat = pd.DataFrame(giuh_cat, columns=['v', 'frequency'])

        giuh_ordinates = ",".join(str(x) for x in np.array(giuh_cat["frequency"]))
        cfe_params.append(f'giuh_ordinates={giuh_ordinates}')
        
        if(runoff_scheme == 'Xinanjiang'):
            cfe_params[1]='surface_partitioning_scheme=Xinanjiang'
            cfe_params.append('a_Xinanjiang_inflection_point_parameter='+str(soil_class_NWM['AXAJ'][soil_id]))
            cfe_params.append('b_Xinanjiang_shape_parameter='+str(soil_class_NWM['BXAJ'][soil_id]))
            cfe_params.append('x_Xinanjiang_shape_parameter='+str(soil_class_NWM['XXAJ'][soil_id]))
            cfe_params.append('urban_decimal_fraction='+str(urban_decimal_fraction))

        # coupled with Soil freeze thaw model
        if(coupled_models == "nom_cfe_smp_sft"):
            cfe_params.append("sft_coupled=true")
            cfe_params.append("ice_content_threshold="+str(ice_content_threshold))

        fname_cfe = cat_name + '_config_cfe.txt'
        cfe_file = os.path.join(cfe_dir, fname_cfe)
        with open(cfe_file, "w") as f:
            f.writelines('\n'.join(cfe_params))

#############################################################################
# The function generates configuration file for TopModel
# @param catids         : array/list of integers contain catchment ids
# @param runoff_schame  : surface runoff schemes - Options = Schaake or Xinanjiang
# @soil_class_NWM       : a dict containing NWM soil characteristics, used here for extracting
#                         Xinanjiang properties for a given soil type
# @gdf_soil             : geodataframe contains soil properties extracted from the hydrofabric
# @param cfe_dir        : output directory (config files are written to this directory)
# @param gpkg_file      : basin geopackage file
# @param coupled_models : option needed to modify CFE config files based on the coupling type
#############################################################################
def write_topmodel_input_files(catids, gdf_soil, topmodel_dir, coupled_models):
    
    # loop over all catchments and write config files
    for catID in catids:
        cat_name = 'cat-'+str(catID) 
        fname = cat_name+'*.txt'
        
        ##################
        topmod = ["0",
                  f'{cat_name}',
                  "./forcing/%s.csv"%cat_name,
                  f'./{topmodel_dir}/subcat_{cat_name}.dat',
                  f'./{topmodel_dir}/params_{cat_name}.dat',
                  f'./{topmodel_dir}/topmod_{cat_name}.out',
                  f'./{topmodel_dir}/hyd_{cat_name}.out'
                  ]

        fname_tm = f'topmod_{cat_name}.run' # + '_config_.run'
        tm_file = os.path.join(topmodel_dir, fname_tm)
        with open(tm_file, "w") as f:
            f.writelines('\n'.join(topmod))

        f.close()
        
        #################
        params = [f'Extracted study basin: {cat_name}',
                  "0.032  5.0  50.  3600.0  3600.0  0.05  0.0000328  0.002  0  1.0  0.02  0.1"
                  ]

        fname_tm = f'params_{cat_name}.dat'
        tm_file = os.path.join(topmodel_dir, fname_tm)
        with open(tm_file, "w") as f:
            f.writelines('\n'.join(params))

        f.close()

        ################
        twi_cat = json.loads(gdf_soil['twi'][cat_name])
        twi_cat = pd.DataFrame(twi_cat, columns=['v', 'frequency'])
        # frequency: distributed area by percentile, v: twi value
                
        twi_cat = twi_cat.sort_values(by=['v'],ascending=False)

        # add width function commulative distribution
        width_f = json.loads(gdf_soil['width_dist'][cat_name])
        df_width_f = pd.DataFrame(width_f, columns=['v', 'frequency'])
        v_cumm = np.cumsum(df_width_f['frequency'])
        
        nclasses_twi = len(twi_cat['frequency'].values)
       
        nclasses_width_function = len(df_width_f['frequency'].values) # width functions (distance to the outlet)
        subcat = ["1 1 1",
                  f'Extracted study basin: {cat_name}',
                  f'{nclasses_twi} 1',
                  'replace_with_twi',
                  f'{nclasses_width_function}',
                  'add_width_function',
                  '$mapfile.dat'
                  ]

        twi_str = ''
        for freq, value in zip(twi_cat['frequency'].values, twi_cat['v'].values):
            twi_str+="{0:.6f}".format(freq) + " " + "{0:.6f}".format(value)+ "\n"
        
        subcat[3] = twi_str.strip()

        # update list/location for the width function
        widthf_str = ''
        for freq, value in zip(v_cumm.values, df_width_f['v'].values):
            widthf_str+="{0:.6f}".format(freq) + " " + "{0:.6f}".format(value)+ " "

        subcat[5] = widthf_str.strip()
        
        fname_tm = f'subcat_{cat_name}.dat'
        tm_file = os.path.join(topmodel_dir, fname_tm)
        with open(tm_file, "w") as f:
            f.writelines('\n'.join(subcat))

        f.close()
                       
#############################################################################
# The function generates configuration file for soil freeze thaw (SFT) model
# @param catids         : array/list of integers contain catchment ids
# @param runoff_schame  : surface runoff schemes - Options = Schaake or Xinanjiang
# @param forcing_dir    : forcing data directory containing data for each catchment
#                         used here for model initialization (based on mean annual air temperature)                
# @param soil_class_NWM : a dict containing NWM soil characteristics, used here for extracting
# @param gdf_soil       : geodataframe contains soil properties extracted from the hydrofabric
#                         Quartz properties for a given soil type
# @param sft_dir        : output directory (config files are written to this directory)
#############################################################################
def write_sft_input_files(catids, runoff_scheme, forcing_dir, gdf_soil, soil_class_NWM, sft_dir):

    # runoff scheme
    if (not runoff_scheme in ["Schaake", "Xinanjiang"]):
        sys.exit("Runoff scheme should be: Schaake or Xinanjiang")

    # num cells -- number of cells used for soil column discretization
    #ncells = 4
    #soil_z = "0.1,0.3,1.0,2.0"

    ncells = 19
    soil_z = "0.1,0.15,0.18,0.23,0.29,0.36,0.44,0.55,0.69,0.86,1.07,1.34,1.66,2.07,2.58,3.22,4.01,5.0,6.0"
    
    delimiter = ','

    nsteps_yr = 365 * 24 # number of steps in the first year of the met. data
    
    # loop over all catchments and write config files
    for catID in catids:
        cat_name = 'cat-'+str(catID)
        forcing_file = glob.glob(os.path.join(forcing_dir, cat_name+'*.csv'))[0]
        
        # obtain annual mean surface temperature as proxy for initial soil temperature
        #df_forcing = pd.read_table(forcing_file,  delimiter=',')
        df_forcing = pd.read_csv(forcing_file,  delimiter=',', usecols=['T2D'], nrows=nsteps_yr, index_col=None)
        
        # compute mean annual air temperature for model initialization
        #MAAT = [str(round(df_forcing['T2D'][:nsteps_yr].mean(), 2)),]*ncells
        MAAT = [str(round(df_forcing['T2D'].mean(), 2)),]*ncells
        MAAT = delimiter.join(MAAT)
        
        # get soil type
        soil_id = gdf_soil['ISLTYP'][cat_name]
        
        # sft params set
        sft_params = ['verbosity=none', 'soil_moisture_bmi=1', 'end_time=1.0[d]', 'dt=1.0[h]', 
                      'soil_params.smcmax=' + str(gdf_soil['soil_params.smcmax'][cat_name]) + '[m/m]', 
                      'soil_params.b=' + str(gdf_soil['soil_params.b'][cat_name]) + '[]', 
                      'soil_params.satpsi=' + str(gdf_soil['soil_params.satpsi'][cat_name]) + '[m]', 
                      'soil_params.quartz=' + str(soil_class_NWM['QTZ'][soil_id]) +'[]',
                      'ice_fraction_scheme=' + runoff_scheme,
                      f'soil_z={soil_z}[m]',
                      f'soil_temperature={MAAT}[K]'
                      ]
         
        fname_sft = cat_name + '_config_sft.txt'
        sft_file = os.path.join(sft_dir, fname_sft)
        with open(sft_file, "w") as f:
            f.writelines('\n'.join(sft_params))

#############################################################################
# The function generates configuration file for soil moisture profiles (SMP) model
# @param catids         : array/list of integers contain catchment ids
# @gdf_soil             : geodataframe contains soil properties extracted from the hydrofabric
# @param smp_dir        : output directory (config files are written to this directory)
# @param coupled_models : option needed to modify SMP config files based on the coupling type
#############################################################################
def write_smp_input_files(catids, gdf_soil, smp_dir, coupled_models):

    #soil_z = "0.1,0.3,1.0,2.0"
    soil_z = "0.1,0.15,0.18,0.23,0.29,0.36,0.44,0.55,0.69,0.86,1.07,1.34,1.66,2.07,2.58,3.22,4.01,5.0,6.0"

    # loop over all catchments and write config files
    for catID in catids:
        cat_name = 'cat-'+str(catID) 

        soil_id = gdf_soil['ISLTYP'][cat_name]
                   
        smp_params = ['verbosity=none',
                      'soil_params.smcmax=' + str(gdf_soil['soil_params.smcmax'][cat_name]) + '[m/m]',
                      'soil_params.b=' + str(gdf_soil['soil_params.b'][cat_name]) + '[]',
                      'soil_params.satpsi=' + str(gdf_soil['soil_params.satpsi'][cat_name]) + '[m]',
                      f'soil_z={soil_z}[m]',
                      'soil_moisture_fraction_depth=1.0[m]'
                      ]

        if ("cfe" in coupled_models):
            smp_params += ['soil_storage_model=conceptual', 'soil_storage_depth=2.0']  
        elif ("lasam" in coupled_models):
            smp_params += ['soil_storage_model=layered', 'soil_moisture_profile_option=constant',
                           'soil_depth_layers=2.0', 'water_table_depth=10[m]']
            # note: soil_depth_layers is an array of depths and will be modified in the future for heterogeneous soils 
            # for exmaple, 'soil_depth_layers=0.4,1.75,2.0'
            # SMCMAX is also an array for hetero. soils

        fname_smp = cat_name + '_config_smp.txt'
        smp_file = os.path.join(smp_dir, fname_smp)
        with open(smp_file, "w") as f:
            f.writelines('\n'.join(smp_params))

#############################################################################
# The function generates configuration file for lumped arid/semi-arid model (LASAM)
# @param catids         : array/list of integers contain catchment ids
# @param giuh_dir       : GIUH data directory (pre-computed GIUH distributions for each catchment)
# @param soil_param_file : input file containing soil properties read by LASAM
#                          (characterizes soil for specified soil types)
# @gdf_soil             : geodataframe contains soil properties extracted from the hydrofabric
# @param lasam_dir        : output directory (config files are written to this directory)
# @param coupled_models : option needed to modify SMP config files based on the coupling type
#############################################################################
def write_lasam_input_files(catids, giuh_dir, soil_param_file, gdf_soil, lasam_dir, coupled_models):

    sft_calib = "False" # update later (should be taken as an argument)

    # used when LASAM coupled with SFT
    #soil_z="10,30,100.0,200.0"
    soil_z = "10.0,15.0,18.0,23.0,29.0,36.0,44.0,55.0,69.0,86.0,107.0,134.0,166.0,207.0,258.0,322.0,401.0,500.0,600.0"
    
    # lasam parameters
    lasam_params_base = ['verbosity=none',
                         'soil_params_file=' + soil_param_file,
                         'layer_thickness=200.0[cm]',
                         'initial_psi=2000.0[cm]',
                         'timestep=300[sec]',
                         'endtime=1000[hr]',
                         'forcing_resolution=3600[sec]',
                         'ponded_depth_max=0[cm]',
                         'use_closed_form_G=false',
                         'layer_soil_type=',
                         'wilting_point_psi=15495.0[cm]',
                         'giuh_ordinates='
                         ]

    if("sft" in coupled_models):
        lasam_params_base.append('sft_coupled=true')
        lasam_params_base.append(f'soil_z={soil_z}[cm]')

    if ( ("sft" in coupled_models) and (sft_calib in ["true", "True"]) ):
        lasam_params_base.append('calib_params=true')

    
    soil_type_loc = lasam_params_base.index("layer_soil_type=")
    giuh_loc_id   = lasam_params_base.index("giuh_ordinates=")
    
    # loop over all catchments and write config files
    for catID in catids:
        cat_name = 'cat-'+str(catID) 
        fname = cat_name+'*.txt'

        #giuh_file = glob.glob(os.path.join(giuh_dir, fname))[0]
        #df_giuh = pd.read_table(giuh_file,  delimiter='=', names=["Params","Values"], index_col=0)

        lasam_params = lasam_params_base.copy()
        lasam_params[soil_type_loc] += str(gdf_soil['ISLTYP'][cat_name])

        # add giuh ordinates
        giuh_cat = json.loads(gdf_soil['giuh'][cat_name])
        giuh_cat = pd.DataFrame(giuh_cat, columns=['v', 'frequency'])

        giuh_ordinates = ",".join(str(x) for x in np.array(giuh_cat["frequency"]))
        lasam_params[giuh_loc_id] += giuh_ordinates

        #lasam_params[giuh_loc_id]  += df_giuh.loc['giuh_ordinates'].iloc[0] # for TauDEM-based workflow
        #lasam_params[giuh_loc_id]   +=  '0.3,0.25,0.2,0.15,0.1'

        
        fname_lasam = cat_name + '_config_lasam.txt'
        lasam_file = os.path.join(lasam_dir, fname_lasam)
        with open(lasam_file, "w") as f:
            f.writelines('\n'.join(lasam_params))

#############################################################################
#############################################################################
def create_directory(dir_name):
    if (os.path.exists(dir_name)):
        str_sub="rm -rf "+dir_name
        out=subprocess.call(str_sub,shell=True)
    os.mkdir(dir_name)


#############################################################################
# main function controlling calls to modules writing config files
# @param gpkg_file : hydrofabric geopackage file (.gpkg)
# @param giuh_dir       : GIUH data directory (pre-computed GIUH distributions for each catchment)
# @param forcing_dir : forcing data directory containing data for each catchment
# @param output_dir        : output directory (config files are written to subdirectories under this directory)
# @param ngen_dir        : path to nextgen directory
# @param models_option : models coupling option (pre-defined names; see main.py)
# @param runoff_schame  : surface runoff schemes - Options = Schaake or Xinanjiang (For CFE and SFT)
# @param time          : dictionary containing simulations start and end time
# @param overwite      : boolean (if true, existing output directories are deleted or overwritten)
#############################################################################
def main():

    try:
        parser = argparse.ArgumentParser()
        parser.add_argument("-gpkg", dest="gpkg_file",     type=str, required=True,  help="the gpkg file")
        parser.add_argument("-giuh", dest="giuh_dir",      type=str, required=False,  help="the giuh files directory")
        parser.add_argument("-f",    dest="forcing_dir",   type=str, required=True,  help="the forcing files directory")
        parser.add_argument("-o",    dest="output_dir",    type=str, required=True,  help="the output files directory")
        parser.add_argument("-ngen", dest="ngen_dir",      type=str, required=True,  help="the ngen directory")
        parser.add_argument("-m",    dest="models_option", type=str, required=True,  help="option for models coupling")
        parser.add_argument("-r",    dest="runoff_scheme", type=str, required=False, help="option for runoff scheme", default="Schaake")
        parser.add_argument("-t",    dest="time",          type=json.loads, required=True,  help="simulation start/end time") 
        parser.add_argument("-ow",   dest="overwrite",     type=str, required=False, default=True,help="overwrite old/existing files")
    except:
        parser.print_help()
        sys.exit(0)
    
    args = parser.parse_args()
    
    #if (not os.path.exists(args.giuh_dir)):
    #    str_msg = 'GIUH directory does not exist! %s'%args.giuh_dir
    #    sys.exit(str_msg)

    if (not os.path.exists(args.gpkg_file)):
        str_msg = 'The gpkg file does not exist! %s'%args.gpkg_file
        sys.exit(str_msg)

    # check if the forcing dir is under Inputs directory
    if (not os.path.exists(args.forcing_dir)):
        str_msg = 'The forcing directory does not exist! %s'%args.forcing_dir
        sys.exit(str_msg)

        
    # giuh dir based cat ids
    #catids = glob.glob(os.path.join(args.giuh_dir,'*.txt'))
    #catids = [int(re.findall('[0-9]+',s)[0]) for s in catids]
   
    gdf_soil, catids = read_gpkg_file(args.gpkg_file, args.models_option)
    
    # doing it outside NOM as some of params from this file are also needed by CFE for Xinanjiang runoff scheme
    nom_params = os.path.join(args.ngen_dir,"extern/noah-owp-modular/noah-owp-modular/parameters")
    
    # *************** NOM  ********************
    if "nom" in args.models_option:
        print ("Generating config files for NOM ...")
        nom_dir = os.path.join(args.output_dir,"nom")
        create_directory(nom_dir)
        str_sub ="cp -r "+ nom_params + " %s"%nom_dir
        out=subprocess.call(str_sub,shell=True)
        nom_soil_file = os.path.join(nom_dir,"parameters/SOILPARM.TBL")

        write_nom_input_files(catids, nom_dir, args.forcing_dir, args.gpkg_file, args.time)
    
    # *************** CFE  ********************
    if "cfe" in args.models_option:
        print ("Generating config files for CFE ...")
        cfe_dir = os.path.join(args.output_dir,"cfe")
        create_directory(cfe_dir)

        # read NWM soil class
        nom_soil_file = os.path.join(nom_params,"SOILPARM.TBL")
        soil_class_NWM = get_soil_class_NWM(nom_soil_file)
        
        write_cfe_input_files(catids, args.runoff_scheme, soil_class_NWM, args.giuh_dir,
                              gdf_soil, cfe_dir, args.models_option)

    # *************** TOPMODEL  ********************
    if "topmodel" in args.models_option:
        print ("Generating config files for TopModel ...")
        tm_dir = os.path.join(args.output_dir,"topmodel")
        create_directory(tm_dir)
        
        write_topmodel_input_files(catids, gdf_soil, tm_dir, args.models_option)
        
    # *************** SFT ********************
    if "sft" in args.models_option:
        print ("Generating config files for SFT and SMP ...")
        smp_only_flag = False
        
        sft_dir = os.path.join(args.output_dir,"sft")
        create_directory(sft_dir)

        smp_dir = os.path.join(args.output_dir,"smp")
        create_directory(smp_dir)

        # read NWM soil class
        nom_soil_file = os.path.join(nom_params,"SOILPARM.TBL")
        soil_class_NWM = get_soil_class_NWM(nom_soil_file)
        
        write_sft_input_files(catids, args.runoff_scheme, args.forcing_dir, gdf_soil,
                              soil_class_NWM, sft_dir)

        write_smp_input_files(catids, gdf_soil, smp_dir, args.models_option)
        
    elif ("smp" in args.models_option):
        print ("Generating config files for SMP...")

        smp_dir = os.path.join(args.output_dir,"smp")
        create_directory(smp_dir)

        write_smp_input_files(catids, gdf_soil, smp_dir, args.models_option)
    
    
    if "lasam" in args.models_option:
        print ("Generating config files for LASAM ...")
        lasam_params = os.path.join(args.ngen_dir,"extern/LGAR-C/data/vG_default_params.dat")

        if (not os.path.isfile(lasam_params)):
            lasam_params = os.path.join(args.ngen_dir,"extern/LASAM/data/vG_default_params.dat")

        lasam_dir = os.path.join(args.output_dir,"lasam")
        create_directory(lasam_dir)
    
        str_sub ="cp -r "+ lasam_params + " %s"%lasam_dir
        out=subprocess.call(str_sub,shell=True)

        write_lasam_input_files(catids, args.giuh_dir, os.path.join(lasam_dir, "vG_default_params.dat"),
                                gdf_soil, lasam_dir, args.models_option)
    
    
    ## create uniform forcings
    #forcing_file = os.path.join(args.forcing_dir,"cat-base.csv")
    #write_forcing_files(catids, forcing_file)
    

if __name__ == "__main__":
    main()

    
