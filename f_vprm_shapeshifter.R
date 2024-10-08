# Function: f_vprm_shapeshifter
# Preparing daily, interpolated VPRM input netcdf files for WRF-GHG 3.9.1.1 runs.
# Interpolation is done treating EVI and LSWI as linear functions of time for each grid cell.
#
# Input:    vprm_input_dir     - directory with MODIS indices (EVI+min+max,
#                                LSWI+min+max, vegetation_fractions) - full path
#                                recommended
#           geo_em_input_dir   - directory with geo_em_dXX files. Will be used
#                                to set the global attributs, necessary for appropriate
#                                handling of the land-use settings. Introduced in v1.2
#           output_dir         - where the output will be stored (full path
#                                recommended)
#           current.domain - character string corresponding to the WRF domain
#                            codes, e.g. c("d01", "d02"). "d01" is default.
#           current.year  - year for which the forecasted VPRM data are being
#                           produced, e.g. 2018 (numeric)
#           previous.year - Leave NULL if NOT working with current calendar year
#                           It's purpose is to handle the last year for which EVI_MAX, EVI_MIN, LSWI_MAX and
#                           LSWI_MIN indices are valid. Usually current.year - 1 (numeric)
#           load.precalculated.indices - debug only
#           no_colons - setting filenames with the no-colon date format, as possible with WRF
#                       namelist option 'nocolons' in the &time_control section.
#                       Default: FALSE
#           add.kaplan.model.input - If input for Kaplan model is available in the input
#                                    dir, it's enabled by this flag
#           kaplan_input_dir - [optional] provide directory with input for Kaplan
#                              will assume vprm_input_dir if not specified
#                              Files need to be located in the input dir and named [units in brackets]:
#                              cpool_fast_LPJ_cdo_dXX.nc [kgC/m2/month],
#                              wetland_kaplan_dXX.nc [fraction (0-1 range)],
#                              t_annual_dXX.nc [K]
# Output:   FILES ONLY    - A year-worth of netcdf files with daily-averaged EVI and
#                           LSWI indices. Will be stored in output_dir


# Debugging:
# On Mistral, one can use these settings for testing purposes:
# vprm_input_dir   <- "/work/mj0143/b301033/Data/CoMet_input/Emissions/VPRM_input/MODIS_indices/JFS.Reanalysis_v2",
#
# geo_em_input_dir <- "/work/mj0143/b301033/Data/CoMet_input/Domains/reanalysis_v2"
# output_dir     <- "/work/mj0143/b301033/Projects/WRF_Tools/vprm_shapeshifter/results"
# requested.domains <- "d01"
# current.year   <- 2018
# previous.year  <- NULL
# load.precalculated.indices <- F
# add.kaplan.model.input     <- F

f_vprm_shapeshifter <- function( vprm_input_dir,
                                 geo_em_input_dir,
                                 output_dir,
                                 requested.domains           = "d01",
                                 current.year                = 2018,
                                 previous.year               = NULL,
				 no_colons                   = FALSE,
                                 load.precalculated.indices  = F,
                                 add.kaplan.model.input      = F,
                                 kaplan_input_dir            = NULL,
                                 kaplan_cpool_file_pattern    = "cpool_fast_LPJ_cdo_dXX.nc", #e.g. cpool_fast_LPJ_cdo_d01.nc
                                 kaplan_wetland_file_pattern  = "wetland_kaplan_dXX.nc",
                                 kaplan_t_annual_file_pattern = "t_annual_dXX.nc" ){
  
  # Previous year paramter is only needed when dealing with a year for which
  # only partial MODIS data is available. The code will then ALSO use previous
  # year's values for MODIS indices (EVI, LSWI) to estimate EVI_MAX, EVI_MIN,
  # LSWI_MAX and LSWI_MIN. Otherwise it is assumed that all data come from
  # a single year.
  if( is.null( previous.year) ){
    previous.year <- current.year
  }
  
  library(ncdf4)
  library(lubridate)
  library(tibble)
  
  cat("\n===================================================================", sep = "")
  cat("\n===================      VPRM Shapeshifter    =====================", sep = "")
  cat("\n===================         for WRF-Chem      =====================", sep = "")
  cat("\n===================           v. 1.5.1        =====================", sep = "")
  cat("\n===================         MPI-BGC 2024      =====================", sep = "")
  cat("\n===================================================================", sep = "")
  cat("\n")
  cat("\nM.Galkowski, MPI-BGC Jena, 2018", sep = "")
  cat("\n")
  cat("\nInput directory:\n")
  cat("  ", vprm_input_dir, "\n", sep = "" )
  cat("\nDirectory with geo_em_dXX files:\n")
  cat("  ", geo_em_input_dir, "\n", sep = "" )
  cat("\nOutput directory:\n")
  cat("  ", output_dir, sep = "" )
  if( !dir.exists( output_dir ) ) {
    cat(" (created)\n" )
    dir.create( output_dir, recursive = TRUE )
  } else {
    cat("\n")
  }
  
  cat("\nDomains requested:\n" )
  cat("  ", paste( requested.domains, collapse = " " ), "\n", sep = "")
  cat("\n")
  
  
  for( current.domain in requested.domains ){
    
    cat("\n-------------------------------------------------------------------", sep = "")
    cat("\nStarting preprocessing of domain ", current.domain, "\n", sep = "")
    cat("\n")
    
    
    
    
    
    # Get the domain configuration from geogrid file ===============================
    cat("Reading domain parameters from: \n")
    path_to_geogrid_file <- file.path( geo_em_input_dir, paste0( "geo_em.", current.domain, ".nc" ) )
    cat("  ", path_to_geogrid_file, "\n")
    
    geo_em_file <- nc_open( path_to_geogrid_file )
    
    lon      <- ncvar_get( geo_em_file, varid = "XLONG_M")
    lat      <- ncvar_get( geo_em_file, varid = "XLAT_M")
    
    n_x <- ncatt_get( geo_em_file, varid = 0, attname = "WEST-EAST_PATCH_END_UNSTAG")$value
    n_y <- ncatt_get( geo_em_file, varid = 0, attname = "SOUTH-NORTH_PATCH_END_UNSTAG")$value
    
    # Get necessary attributes
    geo_em_attributes_floats <- tibble(
      attname = c( "CEN_LAT", "CEN_LON", "TRUELAT1", "TRUELAT2", "MOAD_CEN_LAT", "STAND_LON",
                   "POLE_LAT", "POLE_LON" ),
      value = NA,
      precision = "float"
    )
    for( i in 1:nrow(geo_em_attributes_floats) ){
      geo_em_attributes_floats$value[i] <- ncatt_get( geo_em_file, varid = 0, attname = geo_em_attributes_floats$attname[i] )$value
    }
    
    geo_em_attributes_ints <- tibble(
      attname = c( "MAP_PROJ", "NUM_LAND_CAT", "ISWATER", "ISLAKE", "ISICE", "ISURBAN", "ISOILWATER" ),
      value = NA,
      precision = "int"
    )
    for( i in 1:nrow(geo_em_attributes_ints) ){
      geo_em_attributes_ints$value[i] <- ncatt_get( geo_em_file, varid = 0, attname = geo_em_attributes_ints$attname[i] )$value
    }
    
    geo_em_attribute_MMINLU <- ncatt_get( geo_em_file, varid = 0, attname = "MMINLU" )$value
    
    # DEPRECATED
    # geo_em_atrributes <- tibble(
    #   attname = c("CEN_LAT", "CEN_LON", "TRUELAT1", "TRUELAT2", "MOAD_CEN_LAT", "STAND_LON",
    #               "POLE_LAT", "POLE_LON", "MAP_PROJ",
    #               "MMINLU", "NUM_LAND_CAT", "ISWATER", "ISLAKE", "ISICE", "ISURBAN", "ISOILWATER" ),
    #   # value = as.character("")
    #   value = NA,
    #   precision = c("float", "float", "float", "float", "float", "float", 
    #                 "float", "float", "int",
    #                 "text", "int", "int", "int", "int", "int", "int")
    # )
    # for( i in 1:nrow(geo_em_atrributes) ){
    #   geo_em_atrributes$value[i] <- ncatt_get( geo_em_file, varid = 0, attname = geo_em_atrributes$attname[i] )$value
    # }
    
    nc_close( geo_em_file )
    
    
    
    
    
    # Values of LSWI and EVI will be read from CURRENT YEAR dataset and then
    # interpolated
    lswi.nc <- nc_open( file.path( vprm_input_dir, paste0( "VPRM_input_LSWI_", current.domain, "_", current.year , ".nc" ) ) )
    evi.nc  <- nc_open( file.path( vprm_input_dir, paste0( "VPRM_input_EVI_", current.domain, "_", current.year , ".nc" ) ) )
    
    # Values of min and max indices as well as vegetation fraction will be read from
    # PREVIOUS YEAR datasets
    
    evi.max.nc  <- nc_open( file.path( vprm_input_dir, paste0( "VPRM_input_EVI_MAX_" , current.domain, "_", previous.year , ".nc" ) ) )
    evi.min.nc  <- nc_open( file.path( vprm_input_dir, paste0( "VPRM_input_EVI_MIN_" , current.domain, "_", previous.year , ".nc" ) ) )
    lswi.max.nc <- nc_open( file.path( vprm_input_dir, paste0( "VPRM_input_LSWI_MAX_", current.domain, "_", previous.year , ".nc" ) ) )
    lswi.min.nc <- nc_open( file.path( vprm_input_dir, paste0( "VPRM_input_LSWI_MIN_", current.domain, "_", previous.year , ".nc" ) ) )
    veg_fra.nc  <- nc_open( file.path( vprm_input_dir, paste0( "VPRM_input_VEG_FRA_" , current.domain, "_", previous.year , ".nc" ) ) )
    
    # Added reading of the files necessary for CH4 wetland emissions
    if( add.kaplan.model.input ){
      
      if( is.null( kaplan_input_dir ) ){
        cat("You did not provide separate directory for Kaplan input. We will assume it is located in vprm_input_dir, puny human.\n")
        kaplan_input_dir <- vprm_input_dir
      }
      
      kap.cpool.file  <- sub( "dXX", current.domain, kaplan_cpool_file_pattern )
      kap.wetmap.file <- sub( "dXX", current.domain, kaplan_wetland_file_pattern )
      kap.t.ann.file  <- sub( "dXX", current.domain, kaplan_t_annual_file_pattern )
      
      kap.cpool.nc  <- nc_open( file.path( kaplan_input_dir, kap.cpool.file ) )
      kap.wetmap.nc <- nc_open( file.path( kaplan_input_dir, kap.wetmap.file ) )
      kap.t.ann.nc  <- nc_open( file.path( kaplan_input_dir, kap.t.ann.file ) )
    }
    
    
    # Reading the files into the environment. ====================================
    cat("Reading variables from netcdf files... ", sep = "")
    # 1D
    time <- ncvar_get( lswi.nc, varid = "time" )
    modis.dates <- ISOdate( current.year, 01, 01, 00, tz = "UTC" ) + (time-1)*86400
    
    # 2D ( X x Y )
    lon <- ncvar_get( lswi.nc, varid = "lon")
    lat <- ncvar_get( lswi.nc, varid = "lat")
    
    if( add.kaplan.model.input ){
      kaplan.cpool  <- ncvar_get( kap.cpool.nc,  varid = "CPOOL")  # Expecting units of kgC/m2/month
      kaplan.wetmap <- ncvar_get( kap.wetmap.nc, varid = "WETMAP") # Expecting fraction (0-1 range)!
      kaplan.t.ann  <- ncvar_get( kap.t.ann.nc,  varid = "T_ANN")  # Expecting K
      
      # Dummy fields ( development only, to be erased ):
      # kaplan.cpool   <- rep( 1, length(lon) )
      # kaplan.wetmap  <- rep( 1, length(lon) )
      # kaplan.t.ann   <- rep( 293.15, length(lon) )
    }
    
    # 3D( X x Y x vegetation_class )
    evi.max  <- ncvar_get( evi.max.nc,  varid = "evi_max")
    evi.min  <- ncvar_get( evi.min.nc,  varid = "evi_min")
    lswi.max <- ncvar_get( lswi.max.nc, varid = "lswi_max")
    lswi.min <- ncvar_get( lswi.min.nc, varid = "lswi_min")
    veg.fra  <- ncvar_get( veg_fra.nc,  varid = "vegetation_fraction_map")
    
    # 4D ( X x Y x time x vegetation_class )
    lswi <- ncvar_get( lswi.nc, varid = "lswi" )
    evi  <- ncvar_get( evi.nc,  varid = "evi" )
    cat( "Done.\n", sep = "")
    
    # Closing netcdf connections:
    nc_close( lswi.nc )
    nc_close( evi.nc )
    nc_close( evi.max.nc )
    nc_close( evi.min.nc )
    nc_close( lswi.max.nc )
    nc_close( lswi.min.nc )
    nc_close( veg_fra.nc )
    
    if( add.kaplan.model.input ){
      nc_close( kap.cpool.nc )
      nc_close( kap.wetmap.nc )
      nc_close( kap.t.ann.nc )
    }
    
    
    # Dimensions:
    nx           <- dim(lon)[1]
    ny           <- dim(lon)[2]
    nvegclass    <- dim(lswi)[4]
    
    # TODO: Fix the output dates so that any range can be requested
    # Prepare the output dates:
    #out.dates <- seq( modis.dates[1], rev(modis.dates)[1], by = "day" )
    out.dates <- seq( from = ISOdate( current.year, 01, 01, 00, tz = "UTC"),
                      to   = ISOdate( current.year, 12, 31, 00, tz = "UTC"),
                      by   = "day" )
    
    # This conditional block stores the daily lswi and evi intepolated values as
    # an R object in the temporary directory (hard-set) for later usage. Shortens
    # the run in debugging.
    if( load.precalculated.indices ){
      
      precalculated.indeces.filepath <- paste0( "/work/mj0143/b301033/Projects/CoMet/R/data/daily.evi.and.lswi.", current.domain, ".Rdata" )
      cat( "Precalculated lswi and evi indices were requested (load.precalculated.indices = T).\n" )
      cat( "Please wait, loading from:\n")
      cat("  ", precalculated.indeces.filepath, "\n", sep = "")
      load( precalculated.indeces.filepath )
      cat("Loading completed.\n", sep = "")
      
    } else {
      
      cat( "Requested new interpolation of EVI and LSWI (load.precalculated.indices == F).\n" )
      # Prepare the dummy array that will store intepolated values
      daily.lswi <- array( dim = c( dim( lswi )[1:2], length( out.dates ), nvegclass ) )
      daily.evi  <- array( dim = c( dim( evi )[1:2], length( out.dates ), nvegclass ) )
      cat("Interpolating EVI and LSWI to daily values... \n", sep = "")
      
      # A simple wrapper for apply function to accept the approx function
      # Needed for the apply construct to accept the arguments provided. and
      # to provide output in appropriate format
      f_approx_2d <- function( y, dates.in, dates.out ){
        
        approx( x = dates.in,
                y,
                rule = c(1,2), # rule = c(1,2) returns NA for dates before min(modis.dates), and the nearest value for dates after max(modis.dates)
                xout = dates.out)$y #Approx gives x and y - select y only
        
      }
      
      for( veg.class.idx in 1:nvegclass ){
        
        cat( "  Calculating ncdf vegetation class #", veg.class.idx, "/", n, "\n", sep = "" )
        
        cat( "   Calculating LSWI...\n", sep = "" )
        interpolated.field <- apply( lswi[,,,veg.class.idx], c(1,2), f_approx_2d, dates.in = modis.dates, dates.out = out.dates )
        daily.lswi[,,, veg.class.idx] <- aperm( interpolated.field, c(2,3,1) )
        
        cat( "   Calculating EVI...\n", sep = "" )
        interpolated.field <- apply( evi[,,,veg.class.idx], c(1,2), f_approx_2d, dates.in = modis.dates, dates.out = out.dates )
        daily.evi[,,, veg.class.idx] <- aperm( interpolated.field, c(2,3,1) )
        
      }
      rm( interpolated.field )
      
      # Saving of the interpolated indices as a file to save time
      # Use only when debugging.
      # cat("\nSaving interpolated EVI and LSWI indices. It will take some time, the file can be above 500MB...", sep = "")
      # precalculated.indeces.filepath <- paste0( "/work/mj0143/b301033/Projects/CoMet/R/data/daily.evi.and.lswi.", current.domain, ".Rdata" )
      # save( list = c( "daily.lswi", "daily.evi" ),
      #       file = precalculated.indeces.filepath )
      # cat(" Done.\n", sep = "")
    } # End of if/else with load.precalculated.indices=T/F condition
    
    
    cat("===================================================================\n", sep = "")
    cat("Setting structure of the output netcdf files... ", sep = "")
    
    # Define DIMENSIONS:
    ncdim.Time        <- ncdim_def( name = "Time",        units = "", vals  = 1:1,          create_dimvar = F, unlim = F )
    ncdim.DateStrLen  <- ncdim_def( name = "DateStrLen",  units = "", vals  = 1:19,         create_dimvar = F )
    ncdim.west_east   <- ncdim_def( name = "west_east",   units = "", vals  = 1:nx,         create_dimvar = F )
    ncdim.south_north <- ncdim_def( name = "south_north", units = "", vals  = 1:ny,         create_dimvar = F )
    ncdim.veg_type    <- ncdim_def( name = "vprm_vgcls",  units = "", vals  = 1:nvegclass,  create_dimvar = F )
    ncdim.zdim        <- ncdim_def( name = "zdim",        units = "", vals  = 1:1,          create_dimvar = F )
    
    # DEFINE VARIABLES
    ncvar.times    <- ncvar_def( name = "Times",         units = "",        dim = list( ncdim.DateStrLen, ncdim.Time        ), prec = "char") 
    ncvar.xlong    <- ncvar_def( name = "XLONG",         units = "degrees", dim = list( ncdim.west_east,  ncdim.south_north ), longname = "longitude" )
    ncvar.xlat     <- ncvar_def( name = "XLAT",          units = "degrees", dim = list( ncdim.west_east,  ncdim.south_north ), longname = "latitude" )
    ncvar.evi_min  <- ncvar_def( name = "EVI_MIN",       units = "",        dim = list( ncdim.west_east,  ncdim.south_north,  ncdim.veg_type,  ncdim.Time  ), longname = "minimum annual EVI index value" )
    ncvar.evi_max  <- ncvar_def( name = "EVI_MAX",       units = "",        dim = list( ncdim.west_east,  ncdim.south_north,  ncdim.veg_type,  ncdim.Time  ), longname = "maximum annual EVI index value" )
    ncvar.evi      <- ncvar_def( name = "EVI",           units = "",        dim = list( ncdim.west_east,  ncdim.south_north,  ncdim.veg_type,  ncdim.Time  ), longname = "EVI index value" )
    ncvar.lswi_min <- ncvar_def( name = "LSWI_MIN",      units = "",        dim = list( ncdim.west_east,  ncdim.south_north,  ncdim.veg_type,  ncdim.Time  ), longname = "minimum annual LSWI index value" )
    ncvar.lswi_max <- ncvar_def( name = "LSWI_MAX",      units = "",        dim = list( ncdim.west_east,  ncdim.south_north,  ncdim.veg_type,  ncdim.Time  ), longname = "maximum annual LSWI index value" )
    ncvar.lswi     <- ncvar_def( name = "LSWI",          units = "",        dim = list( ncdim.west_east,  ncdim.south_north,  ncdim.veg_type,  ncdim.Time  ), longname = "LSWI index value" )
    ncvar.vegfra   <- ncvar_def( name = "VEGFRA_VPRM",   units = "",        dim = list( ncdim.west_east,  ncdim.south_north,  ncdim.veg_type,  ncdim.Time  ), longname = "VPRM vegetation fraction" )
    # DEV:
    # Additional dummy variables for Kaplan model that WRF might want.
    # Comment out if not needed
    if( add.kaplan.model.input ){
      ncvar.kaplan.cpool   <- ncvar_def( name = "CPOOL",   units = "gC/m^2",  dim = list( ncdim.west_east, ncdim.south_north, ncdim.zdim, ncdim.Time ), longname = "LPJ Carbon pool" )
      ncvar.kaplan.wetmap  <- ncvar_def( name = "WETMAP",  units = "",        dim = list( ncdim.west_east, ncdim.south_north, ncdim.zdim, ncdim.Time ), longname = "Kaplan potential wetland map" )
      ncvar.kaplan.t.ann   <- ncvar_def( name = "T_ANN",   units = "K",       dim = list( ncdim.west_east, ncdim.south_north, ncdim.zdim, ncdim.Time ), longname = "mean annual temperature" )
    }
    
    cat("Done \n", sep = "")
    

    # Writing output: ============================================================
    cat("Writing VPRM daily interpolated files...\n", sep = "")
    
    n.times <- length(out.dates)
    
    # New in version 1.5: choose filename date format
    filename_date_format <- if( no_colons ){ "%Y-%m-%d_%H_%M_%S" } else { "%Y-%m-%d_%H:%M:%S" }

    for( time.idx in 1:n.times ){
      
      current.date.code <- format( out.dates[time.idx], format = "%Y-%m-%d_%H:%M:%S" )
      current.date.code.for.filename <- format( out.dates[time.idx], format = filename_date_format )
      current.filename  <- paste0( "vprm_input_", current.domain, "_", current.date.code.for.filename, ".nc" )
      
      cat( "\r #", time.idx, "/", n.times, ": ", current.filename, "", sep = "" )
      flush.console()
      
      # CREATE A LIST OF ALL VARIABLES
      list.of.variables <- list( ncvar.times,
                                 ncvar.xlong, ncvar.xlat,
                                 ncvar.evi_min, ncvar.evi_max, ncvar.evi,
                                 ncvar.lswi_min, ncvar.lswi_max, ncvar.lswi,
                                 ncvar.vegfra )
      
      if( add.kaplan.model.input ){
        list.of.variables <- c( list.of.variables,
                                list( ncvar.kaplan.cpool, ncvar.kaplan.wetmap, ncvar.kaplan.t.ann) )
      }
      
      # CREATE NCDF FILE AND ASSIGN VALUES TO VARIABLES
      ncnew <- nc_create( filename = file.path( output_dir, current.filename ),
                          vars     = list.of.variables,
                          force_v4 = TRUE)
      
      ncvar_put( nc = ncnew, varid = ncvar.times,    vals = current.date.code )
      ncvar_put( nc = ncnew, varid = ncvar.xlong,    vals = as.vector( lon ) )
      ncvar_put( nc = ncnew, varid = ncvar.xlat,     vals = as.vector( lat ) )
      ncvar_put( nc = ncnew, varid = ncvar.evi_min,  vals = as.vector( evi.min ) )
      ncvar_put( nc = ncnew, varid = ncvar.evi_max,  vals = as.vector( evi.max ) )
      ncvar_put( nc = ncnew, varid = ncvar.evi,      vals = as.vector( aperm( daily.evi, perm = c(1,2,4,3) )[,,,time.idx] ) )
      ncvar_put( nc = ncnew, varid = ncvar.lswi_min, vals = as.vector( lswi.min ) )
      ncvar_put( nc = ncnew, varid = ncvar.lswi_max, vals = as.vector( lswi.max ) )
      ncvar_put( nc = ncnew, varid = ncvar.lswi,     vals = as.vector( aperm( daily.lswi, perm = c(1,2,4,3) )[,,,time.idx] ) )
      ncvar_put( nc = ncnew, varid = ncvar.vegfra,   vals = as.vector( veg.fra ) )
      
      # DEV: Dummy output for Kaplan model.
      if( add.kaplan.model.input ){
        ncvar_put( nc = ncnew, varid = ncvar.kaplan.cpool,   vals = kaplan.cpool  )
        ncvar_put( nc = ncnew, varid = ncvar.kaplan.wetmap,  vals = kaplan.wetmap )
        ncvar_put( nc = ncnew, varid = ncvar.kaplan.t.ann,   vals = kaplan.t.ann  )
      }
      
      
      
      
      
      # GLOBAL ATTRIBUTES (varid = 0 means global):
      # These are NOT optional. Not having those attributes set can cause issues
      # with the land-use settings during simulation!!!
      # Michal Galkowski, April 2019
      
      # These are NECESSARY, because WRF actually reads those. And not only reads,
      # but landuse information is the one that drives the model.
      # Mike G - see personal logfile entry from March 8, 2019
      ncatt_put( nc = ncnew, varid = 0, attname = "Source", attval = "WRF input file created by vprm_shapeshifter, v1.3 (MPI-BGC Jena 2019)")
      
      
      for( i in 1:nrow(geo_em_attributes_floats) ){
        ncatt_put( nc = ncnew,
                   varid = 0,
                   attname = geo_em_attributes_floats$attname[i],
                   attval  = geo_em_attributes_floats$value[i],
                   prec    = geo_em_attributes_floats$precision[i] )
      } # End of writing attributes from geo_em_attributes_floats
      
      for( i in 1:nrow(geo_em_attributes_ints) ){
        ncatt_put( nc = ncnew,
                   varid = 0,
                   attname = geo_em_attributes_ints$attname[i],
                   attval  = geo_em_attributes_ints$value[i],
                   prec    = geo_em_attributes_ints$precision[i] )
      } # End of writing attributes from geo_em_attributes_ints
      
      # Now, write also additional time-dependent attributes.
      ncatt_put( nc = ncnew, varid = 0, attname = "MMINLU", attval = geo_em_attribute_MMINLU, prec = "text" )
      ncatt_put( nc = ncnew, varid = 0, attname = "GMT", attval = 0, prec = "float")
      ncatt_put( nc = ncnew, varid = 0, attname = "JULYR", attval = year( out.dates[time.idx] ), prec = "int" )
      ncatt_put( nc = ncnew, varid = 0, attname = "JULDAY", attval = yday( out.dates[time.idx] ), prec = "int" )

      
      
      
      
      # VARIABLE ATTRIBUTES:
      ncatt_put( nc = ncnew, varid = ncvar.times, attname = "description", attval = "WRF-format date_time string" )
      
      ncatt_put( nc = ncnew, varid = ncvar.xlong, attname = "FieldType",   attval = 104, prec = "int" )
      ncatt_put( nc = ncnew, varid = ncvar.xlong, attname = "MemoryOrder", attval = "XY " )
      ncatt_put( nc = ncnew, varid = ncvar.xlong, attname = "units",       attval = "degrees longitude" )
      ncatt_put( nc = ncnew, varid = ncvar.xlong, attname = "description", attval = "Longitude on mass grid" )
      ncatt_put( nc = ncnew, varid = ncvar.xlong, attname = "stagger",     attval = "M" )
      
      ncatt_put( nc = ncnew, varid = ncvar.xlat, attname = "FieldType",    attval = 104, prec = "int" )
      ncatt_put( nc = ncnew, varid = ncvar.xlat, attname = "MemoryOrder",  attval = "XY " )
      ncatt_put( nc = ncnew, varid = ncvar.xlat, attname = "units",        attval = "degrees latitude" )
      ncatt_put( nc = ncnew, varid = ncvar.xlat, attname = "description",  attval = "Latitude on mass grid" )
      ncatt_put( nc = ncnew, varid = ncvar.xlat, attname = "stagger",      attval = "M" )
      
      ncatt_put( nc = ncnew, varid = ncvar.evi_min, attname = "FieldType",   attval = 104, prec = "int" )
      ncatt_put( nc = ncnew, varid = ncvar.evi_min, attname = "MemoryOrder", attval = "XYZ" )
      ncatt_put( nc = ncnew, varid = ncvar.evi_min, attname = "units",       attval = "" )
      ncatt_put( nc = ncnew, varid = ncvar.evi_min, attname = "description", attval = "Minimal value of EVI index" )
      ncatt_put( nc = ncnew, varid = ncvar.evi_min, attname = "stagger",     attval = "M" )
      ncatt_put( nc = ncnew, varid = ncvar.evi_min, attname = "coordinates", attval = "XLONG XLAT" )
      
      ncatt_put( nc = ncnew, varid = ncvar.evi_max, attname = "FieldType",   attval = 104, prec = "int" )
      ncatt_put( nc = ncnew, varid = ncvar.evi_max, attname = "MemoryOrder", attval = "XYZ" )
      ncatt_put( nc = ncnew, varid = ncvar.evi_max, attname = "units",       attval = "" )
      ncatt_put( nc = ncnew, varid = ncvar.evi_max, attname = "description", attval = "Maximal value of EVI index" )
      ncatt_put( nc = ncnew, varid = ncvar.evi_max, attname = "stagger",     attval = "M" )
      ncatt_put( nc = ncnew, varid = ncvar.evi_max, attname = "coordinates", attval = "XLONG XLAT" )
      
      ncatt_put( nc = ncnew, varid = ncvar.evi, attname = "FieldType",   attval = 104, prec = "int" )
      ncatt_put( nc = ncnew, varid = ncvar.evi, attname = "MemoryOrder", attval = "XYZ" )
      ncatt_put( nc = ncnew, varid = ncvar.evi, attname = "units",       attval = "" )
      ncatt_put( nc = ncnew, varid = ncvar.evi, attname = "description", attval = "Value of EVI index" )
      ncatt_put( nc = ncnew, varid = ncvar.evi, attname = "stagger",     attval = "M" )
      ncatt_put( nc = ncnew, varid = ncvar.evi, attname = "coordinates", attval = "XLONG XLAT" )
      
      ncatt_put( nc = ncnew, varid = ncvar.lswi_min, attname = "FieldType",   attval = 104, prec = "int" )
      ncatt_put( nc = ncnew, varid = ncvar.lswi_min, attname = "MemoryOrder", attval = "XYZ" )
      ncatt_put( nc = ncnew, varid = ncvar.lswi_min, attname = "units",       attval = "" )
      ncatt_put( nc = ncnew, varid = ncvar.lswi_min, attname = "description", attval = "Minimal value of lswi index" )
      ncatt_put( nc = ncnew, varid = ncvar.lswi_min, attname = "stagger",     attval = "M" )
      ncatt_put( nc = ncnew, varid = ncvar.lswi_min, attname = "coordinates", attval = "XLONG XLAT" )
      
      ncatt_put( nc = ncnew, varid = ncvar.lswi_max, attname = "FieldType",   attval = 104, prec = "int" )
      ncatt_put( nc = ncnew, varid = ncvar.lswi_max, attname = "MemoryOrder", attval = "XYZ" )
      ncatt_put( nc = ncnew, varid = ncvar.lswi_max, attname = "units",       attval = "" )
      ncatt_put( nc = ncnew, varid = ncvar.lswi_max, attname = "description", attval = "Maximal value of lswi index" )
      ncatt_put( nc = ncnew, varid = ncvar.lswi_max, attname = "stagger",     attval = "M" )
      ncatt_put( nc = ncnew, varid = ncvar.lswi_max, attname = "coordinates", attval = "XLONG XLAT" )
      
      ncatt_put( nc = ncnew, varid = ncvar.lswi, attname = "FieldType",   attval = 104, prec = "int" )
      ncatt_put( nc = ncnew, varid = ncvar.lswi, attname = "MemoryOrder", attval = "XYZ" )
      ncatt_put( nc = ncnew, varid = ncvar.lswi, attname = "units",       attval = "" )
      ncatt_put( nc = ncnew, varid = ncvar.lswi, attname = "description", attval = "Value of lswi index" )
      ncatt_put( nc = ncnew, varid = ncvar.lswi, attname = "stagger",     attval = "M" )
      ncatt_put( nc = ncnew, varid = ncvar.lswi, attname = "coordinates", attval = "XLONG XLAT" )
      
      ncatt_put( nc = ncnew, varid = ncvar.vegfra, attname = "FieldType",   attval = 104, prec = "int" )
      ncatt_put( nc = ncnew, varid = ncvar.vegfra, attname = "MemoryOrder", attval = "XYZ" )
      ncatt_put( nc = ncnew, varid = ncvar.vegfra, attname = "units",       attval = "" )
      ncatt_put( nc = ncnew, varid = ncvar.vegfra, attname = "description", attval = "Vegetation fraction for VPRM" )
      ncatt_put( nc = ncnew, varid = ncvar.vegfra, attname = "stagger",     attval = "M" )
      ncatt_put( nc = ncnew, varid = ncvar.vegfra, attname = "coordinates", attval = "XLONG XLAT" )
      
      
      # WRF-GHG Kaplan model input
      if( add.kaplan.model.input ){
        ncatt_put( nc = ncnew, varid = ncvar.kaplan.cpool, attname = "FieldType", attval = 104, prec = "int" )
        ncatt_put( nc = ncnew, varid = ncvar.kaplan.cpool, attname = "MemoryOrder", attval = "XYZ" )
        ncatt_put( nc = ncnew, varid = ncvar.kaplan.cpool, attname = "units", attval = "kgC/m^2" )
        ncatt_put( nc = ncnew, varid = ncvar.kaplan.cpool, attname = "description", attval = "Carbon pool value for Kaplan model" )
        ncatt_put( nc = ncnew, varid = ncvar.kaplan.cpool, attname = "stagger", attval = "M" )
        ncatt_put( nc = ncnew, varid = ncvar.kaplan.cpool, attname = "coordinates", attval = "XLONG XLAT" )
        
        ncatt_put( nc = ncnew, varid = ncvar.kaplan.wetmap, attname = "FieldType", attval = 104, prec = "int" )
        ncatt_put( nc = ncnew, varid = ncvar.kaplan.wetmap, attname = "MemoryOrder", attval = "XYZ" )
        ncatt_put( nc = ncnew, varid = ncvar.kaplan.wetmap, attname = "units", attval = "1 Woolong" )
        ncatt_put( nc = ncnew, varid = ncvar.kaplan.wetmap, attname = "description", attval = "Wetland map for Kaplan model" )
        ncatt_put( nc = ncnew, varid = ncvar.kaplan.wetmap, attname = "stagger", attval = "M" )
        ncatt_put( nc = ncnew, varid = ncvar.kaplan.wetmap, attname = "coordinates", attval = "XLONG XLAT" )
        
        ncatt_put( nc = ncnew, varid = ncvar.kaplan.t.ann, attname = "FieldType", attval = 104, prec = "int" )
        ncatt_put( nc = ncnew, varid = ncvar.kaplan.t.ann, attname = "MemoryOrder", attval = "XYZ" )
        ncatt_put( nc = ncnew, varid = ncvar.kaplan.t.ann, attname = "units", attval = "K" )
        ncatt_put( nc = ncnew, varid = ncvar.kaplan.t.ann, attname = "description", attval = "Annual mean temperature for vegetation classes" )
        ncatt_put( nc = ncnew, varid = ncvar.kaplan.t.ann, attname = "stagger", attval = "M" )
        ncatt_put( nc = ncnew, varid = ncvar.kaplan.t.ann, attname = "coordinates", attval = "XLONG XLAT" )
      }
      
      nc_close( ncnew )
      
      
    } # End of for loop with time.idx
    
    
  } # End of a single domain preprocessing block
  
  cat("\n=========================== THE END ===============================\n\n\n", sep = "")
  
} # End of function
