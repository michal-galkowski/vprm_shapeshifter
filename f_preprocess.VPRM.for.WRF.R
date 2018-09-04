# Function: f_preprocess.VPRM.for.WRF
# Preparing daily, interpolated VPRM input netcdf files for WRF-GHG 3.9.1.1 runs.
# Interpolation is done treating EVI and LSWI as linear functions of time for each grid cell.
#
# Input:    input_dir     - directory with MODIS indices (EVI+min+max,
#                           LSWI+min+max, vegetation_fractions) - full path
#                           recommended
#           output_dir    - where the output will be stored (full path
#                           recommended)
#           current.domain - character string corresponding to the WRF domain
#                           code, e.g. "d01", "d02"...
#           current.year  - year for which the forecasted VPRM data are being
#                           produced, e.g. 2018 (numeric)
#           previous.year - last year for which EVI_MAX. EVI_MIN, LSWI_MAX and
#                           LSWI_MIN indices are valid. E.g. 2017 (numeric)
#
# Output:   FILES ONLY    - Bunch of netcdf files with daily-averaged EVI and
#                           LSWI indices.


# Debugging:
# On Mistral, one can use these settings for testing purposes:
# input_dir      = "/work/mj0143/b301033/Data/CoMet_input/Emissions/VPRM_input/yearly/Operational.Data"
# output_dir     = "/work/mj0143/b301033/Data/CoMet_input/Emissions/VPRM_input/WRF_input"
# current.domain <- "d01"
# current.year   <- 2018
# previous.year  <- 2017

f_preprocess.VPRM.for.WRF <- function( input_dir                  = "/work/mj0143/b301033/Data/CoMet_input/Emissions/VPRM_input/yearly",
                                       output_dir                 = "/work/mj0143/b301033/Data/CoMet_input/Emissions/VPRM_input/WRF_input",
                                       current.domain             = "d01",
                                       current.year               = 2018,
                                       previous.year              = 2017,
                                       load.precalculated.indices = F ){
  
  library(ncdf4)
  
  cat("\n===================================================================", sep = "")
  cat("\n=================== VPRM for WRF-Chem 3.9.1.1.=====================", sep = "")
  cat("\n===================         CoMet 2018        =====================", sep = "")
  cat("\n===================================================================", sep = "")
  cat("\n")
  cat("Input directory:\n")
  cat("  ", input_dir, "\n", sep = "")
  cat("Outputt directory:\n")
  cat("  ", output_dir, "\n", sep = "")
  cat("Selected domain:\n")
  cat("  ", current.domain, "\n\n", sep = "")
  
  # Values of LSWI and EVI will be read from CURRENT YEAR dataset and then
  # interpolated
  lswi.nc <- nc_open( file.path( input_dir, paste0( "VPRM_input_LSWI_", current.domain, "_", current.year , ".nc" ) ) )
  evi.nc  <- nc_open( file.path( input_dir, paste0( "VPRM_input_EVI_", current.domain, "_", current.year , ".nc" ) ) )
  
  # Values of min and max indices as well as vegetation fraction will be read from
  # PREVIOUS YEAR datasets
  
  evi.max.nc  <- nc_open( file.path( input_dir, paste0( "VPRM_input_EVI_MAX_" , current.domain, "_", previous.year , ".nc" ) ) )
  evi.min.nc  <- nc_open( file.path( input_dir, paste0( "VPRM_input_EVI_MIN_" , current.domain, "_", previous.year , ".nc" ) ) )
  lswi.max.nc <- nc_open( file.path( input_dir, paste0( "VPRM_input_LSWI_MAX_", current.domain, "_", previous.year , ".nc" ) ) )
  lswi.min.nc <- nc_open( file.path( input_dir, paste0( "VPRM_input_LSWI_MIN_", current.domain, "_", previous.year , ".nc" ) ) )
  veg_fra.nc  <- nc_open( file.path( input_dir, paste0( "VPRM_input_VEG_FRA_" , current.domain, "_", previous.year , ".nc" ) ) )
  
  # Reading the files into the environment. ====================================
  cat("Reading variables from netcdf files... ", sep = "")
  # 1D
  time <- ncvar_get( lswi.nc, varid = "time" )
  date <- ISOdate( current.year, 01, 01, 00, tz = "UTC" ) + (time-1)*86400
  
  # 2D ( X x Y )
  lon <- ncvar_get( lswi.nc, varid = "lon")
  lat <- ncvar_get( lswi.nc, varid = "lat")
  
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
  
  # Attributes from a met_em file from one of the runs.
  # Copying is faster than writing them from scratch.
  # TODO: Do it anyway, creating a dummy file is suboptimal
  # attributes.nc        <- nc_open( file.path( input_dir, paste0( "attributes_dummy_file_", current.domain, ".nc" ) ) )
  # attributes.values    <- ncatt_get( attributes.nc, varid = 0, attname = NA )
  
  # Dummy fields:
  kaplan.cpool   <- rep( 1, length(veg.fra) )
  kaplan.wetmap  <- rep( 1, length(lon)*8 )
  kaplan.t.ann   <- rep( 293.15, length(lon)*8 )
  
  
  # Closing netcdf connections:
  nc_close( lswi.nc )
  nc_close( evi.nc )
  nc_close( evi.max.nc )
  nc_close( evi.min.nc )
  nc_close( lswi.max.nc )
  nc_close( lswi.min.nc )
  nc_close( veg_fra.nc )
  # nc_close( attributes.nc )
  
  
  # Dimensions:
  nx <- dim(lon)[1]
  ny <- dim(lon)[2]
  
  # Prepare the output dates:
  out.dates <- seq( date[1], rev(date)[1], by = "day" )
  
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
    
    cat( "Requested new interpolation of EVI and LSWI (load.precalculated.indices = F).\n" )
    # Prepare the dummy array that will store intepolated values
    daily.lswi <- array( dim = c( dim( lswi )[1:2], length( out.dates ), dim( lswi )[4] ) )
    daily.evi  <- array( dim = c( dim( evi )[1:2], length( out.dates ), dim( evi )[4] ) )
    n <- dim(daily.lswi)[1]
    
    cat("Interpolating EVI and LSWI to daily values... \n", sep = "")
    
    # Super ugly for-loop structure:
    for( i in 1:dim(daily.lswi)[1]){
      cat( "\r  Calculating ncdf longitude band #", i, "/", n, "", sep = "" )
      for( j in 1:dim(daily.lswi)[2]){
        for( veg.class.idx in 1:dim(daily.lswi)[4]){
          
          x <- date
          y <- lswi[ i, j, ,veg.class.idx ]
          daily.lswi[ i, j, , veg.class.idx ] <- approx( x, y, xout = out.dates )$y
          
          x <- date
          y <- evi[ i, j, ,veg.class.idx ]
          daily.evi[ i, j, , veg.class.idx ] <- approx( x, y, xout = out.dates )$y
          
        } # End of for loop with veg.class.idx index.
      } # End of for loop with j index.
    } # End of for loop with i index.
    
    cat("\nSaving interpolated EVI and LSWI indices. It will take some time, the file can be above 500MB...", sep = "")
    precalculated.indeces.filepath <- paste0( "/work/mj0143/b301033/Projects/CoMet/R/data/daily.evi.and.lswi.", current.domain, ".Rdata" )
    save( list = c( "daily.lswi", "daily.evi" ),
          file = precalculated.indeces.filepath )
    cat(" Done.\n", sep = "")
  } # End of if/else with load.precalculated.indices=T/F condition
  
  
  cat("===================================================================\n", sep = "")
  cat("Setting structure of the output netcdf files... ", sep = "")
  
  # Define DIMENSIONS:
  ncdim.Time        <- ncdim_def( name = "Time",        units = "", vals  = 1:1,  create_dimvar = F, unlim = F )
  ncdim.DateStrLen  <- ncdim_def( name = "DateStrLen",  units = "", vals  = 1:19, create_dimvar = F )
  ncdim.west_east   <- ncdim_def( name = "west_east",   units = "", vals  = 1:nx, create_dimvar = F )
  ncdim.south_north <- ncdim_def( name = "south_north", units = "", vals  = 1:ny, create_dimvar = F )
  ncdim.veg_type    <- ncdim_def( name = "vprm_vgcls",  units = "", vals  = 1:8,  create_dimvar = F )
  
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
  ncvar.kaplan.cpool   <- ncvar_def( name = "CPOOL",   units = "gC/m^2",  dim = list( ncdim.west_east, ncdim.south_north, ncdim.veg_type, ncdim.Time ), longname = "LPJ Carbon pool" )
  ncvar.kaplan.wetmap  <- ncvar_def( name = "WETMAP",  units = "",        dim = list( ncdim.west_east, ncdim.south_north, ncdim.veg_type, ncdim.Time ), longname = "Kaplan potential wetland map" )
  ncvar.kaplan.t.ann   <- ncvar_def( name = "T_ANN",   units = "K",       dim = list( ncdim.west_east, ncdim.south_north, ncdim.veg_type, ncdim.Time ), longname = "mean annual temperature" )
  
  
  # Writing output: ============================================================
  cat("Done \n", sep = "")
  cat("Writing VPRM daily interpolated files...\n", sep = "")
  
  n.times <- length(out.dates)
  
  for( time.idx in 1:n.times ){
    
    current.date.code <- format( out.dates[time.idx], format = "%Y-%m-%d_%H:%M:%S" )
    current.filename  <- paste0( "vprm_input_", current.domain, "_", current.date.code )
    
    cat( "\r #", time.idx, "/", n.times, ": ", current.filename, "", sep = "" )
    flush.console()
    
    # CREATE NCDF FILE AND ASSIGN VALUES TO VARIABLES
    ncnew <- nc_create( filename = file.path( output_dir, current.filename ),
                        vars     = list( ncvar.times,
                                         ncvar.xlong, ncvar.xlat,
                                         ncvar.evi_min, ncvar.evi_max, ncvar.evi,
                                         ncvar.lswi_min, ncvar.lswi_max, ncvar.lswi,
                                         ncvar.vegfra,
                                         ncvar.kaplan.cpool, ncvar.kaplan.wetmap, ncvar.kaplan.t.ann),
                        force_v4 = F)
    
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
    ncvar_put( nc = ncnew, varid = ncvar.kaplan.cpool,   vals = kaplan.cpool  )
    ncvar_put( nc = ncnew, varid = ncvar.kaplan.wetmap,  vals = kaplan.wetmap )
    ncvar_put( nc = ncnew, varid = ncvar.kaplan.t.ann,   vals = kaplan.t.ann  )
    
    
    
    # GLOBAL ATTRIBUTES (varid = 0 means global):
    # Right now copying from wrfinput/met_em files.
    # Optimally change that to include only what is necessary
    # for( i in 1:length(attributes.values) ){
    #   ncatt_put( nc = ncnew, varid = 0, attname = names( attributes.values[i] ), attval = attributes.values[[i]] )
    # }
    
    # ncatt_put( nc = ncnew, varid = 0, attname = "CEN_LAT" )
    # ncatt_put( nc = ncnew, varid = 0, attname = "CEN_LON" )
    # ncatt_put( nc = ncnew, varid = 0, attname = "TRUELAT1" )
    # ncatt_put( nc = ncnew, varid = 0, attname = "TRUELAT2" )
    # ncatt_put( nc = ncnew, varid = 0, attname = "MOAD_CEN_LAT" )
    # ncatt_put( nc = ncnew, varid = 0, attname = "STAND_LON" )
    # ncatt_put( nc = ncnew, varid = 0, attname = "POLE_LAT" )
    # ncatt_put( nc = ncnew, varid = 0, attname = "POLE_LON" )
    # ncatt_put( nc = ncnew, varid = 0, attname = "GMT" )
    # ncatt_put( nc = ncnew, varid = 0, attname = "JULYR" )
    # ncatt_put( nc = ncnew, varid = 0, attname = "JULDAY" )
    # ncatt_put( nc = ncnew, varid = 0, attname = "MAP_PROJ" )
    # ncatt_put( nc = ncnew, varid = 0, attname = "MMINLU" )
    # ncatt_put( nc = ncnew, varid = 0, attname = "ISWATER" )
    # ncatt_put( nc = ncnew, varid = 0, attname = "ISLAKE" )
    # ncatt_put( nc = ncnew, varid = 0, attname = "ISICE" )
    # ncatt_put( nc = ncnew, varid = 0, attname = "ISURBAN" )
    # ncatt_put( nc = ncnew, varid = 0, attname = "ISOILWATER" )
    
    
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
    
    
    # DEV: Dummy output for Kaplan model
    ncatt_put( nc = ncnew, varid = ncvar.kaplan.cpool, attname = "FieldType", attval = 104, prec = "int" )
    ncatt_put( nc = ncnew, varid = ncvar.kaplan.cpool, attname = "MemoryOrder", attval = "XYZ" )
    ncatt_put( nc = ncnew, varid = ncvar.kaplan.cpool, attname = "units", attval = "gC/m^2" )
    ncatt_put( nc = ncnew, varid = ncvar.kaplan.cpool, attname = "description", attval = "Carbon pool value for Kaplan model" )
    ncatt_put( nc = ncnew, varid = ncvar.kaplan.cpool, attname = "stagger", attval = "M" )
    ncatt_put( nc = ncnew, varid = ncvar.kaplan.cpool, attname = "coordinates", attval = "XLONG XLAT" )
    
    ncatt_put( nc = ncnew, varid = ncvar.kaplan.wetmap, attname = "FieldType", attval = 104, prec = "int" )
    ncatt_put( nc = ncnew, varid = ncvar.kaplan.wetmap, attname = "MemoryOrder", attval = "XYZ" )
    ncatt_put( nc = ncnew, varid = ncvar.kaplan.wetmap, attname = "units", attval = "something something dark side" )
    ncatt_put( nc = ncnew, varid = ncvar.kaplan.wetmap, attname = "description", attval = "Wetland map for Kaplan model" )
    ncatt_put( nc = ncnew, varid = ncvar.kaplan.wetmap, attname = "stagger", attval = "M" )
    ncatt_put( nc = ncnew, varid = ncvar.kaplan.wetmap, attname = "coordinates", attval = "XLONG XLAT" )
    
    ncatt_put( nc = ncnew, varid = ncvar.kaplan.t.ann, attname = "FieldType", attval = 104, prec = "int" )
    ncatt_put( nc = ncnew, varid = ncvar.kaplan.t.ann, attname = "MemoryOrder", attval = "XYZ" )
    ncatt_put( nc = ncnew, varid = ncvar.kaplan.t.ann, attname = "units", attval = "K" )
    ncatt_put( nc = ncnew, varid = ncvar.kaplan.t.ann, attname = "description", attval = "Annual mean temperature for vegetation classes" )
    ncatt_put( nc = ncnew, varid = ncvar.kaplan.t.ann, attname = "stagger", attval = "M" )
    ncatt_put( nc = ncnew, varid = ncvar.kaplan.t.ann, attname = "coordinates", attval = "XLONG XLAT" )
    
    nc_close( ncnew )
    
    
  } # End of for loop with time.idx
  
  cat("\n=========================== THE END ===============================\n\n\n", sep = "")
  return( 0 )
} # End of function
