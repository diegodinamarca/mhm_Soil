{
  # HEADER --------------------------------------------
  #
  # Author: Diego Dinamarca
  # Email:  ddinamarcamuller@gmail.com
  # 
  # Date:
  #
  # Script Name:
  #
  # Script Description:
  #
  # Notes:
  #
  #
  # INSTALL PACKAGES & LOAD LIBRARIES -----------------
  cat("INSTALLING PACKAGES & LOADING LIBRARIES... \n\n", sep = "")
  packages <- c("tidyverse",
                "terra", 
                "sf",
                "here",
                "gridExtra",
                "raster",
                "rasterVis",
                "ggpointdensity",
                "viridis",
                "scico",
                "ncdf4",
                "chron",
                "stars") # list of packages to load
  n_packages <- length(packages) # count how many packages are required
  
  new.pkg <- packages[!(packages %in% installed.packages())] # determine which packages aren't installed
  
  # install missing packages
  if(length(new.pkg)){
    install.packages(new.pkg)
  }
  
  # load all requried libraries
  for(n in 1:n_packages){
    cat("Loading Library #", n, " of ", n_packages, "... Currently Loading: ", packages[n], "\n", sep = "")
    lib_load <- paste("library(\"",packages[n],"\")", sep = "") # create string of text for loading each library
    eval(parse(text = lib_load)) # evaluate the string to load the library
  }
  # SET WORKING DIRECTORY -----------------------------
  cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
  wd <- here::here()
  setwd(wd)
  cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")
  
  # SET OPTIONS ---------------------------------------
  cat("SETTING OPTIONS... \n\n", sep = "")
  options(scipen = 999) # turns off scientific notation
  options(encoding = "UTF-8") # sets string encoding to UTF-8 instead of ANSI
  
  
  # CONFLICTS ---------------------------------------------------------------
  conflicted::conflict_prefer("select", "dplyr")
  conflicted::conflict_prefer("filter", "dplyr")
  conflicted::conflict_prefer("extract", "terra")
  
  # LOAD FUNCTIONS ------------------------------------
}

m = "CLSM"
# restart grids directory
models = c("CLSM","SG","FAO")
variables = c("L1_soilMoistFC", "L1_wiltingPoint")
for (k in 1:length(models)) {
  # k=1
  
  m = models[k]
  n = length(variables)
  dir = here("BASIN_DATA/out/7339001/calibration_runs",str_c(m,"_calibrated"))
  files = list.files(dir, recursive = TRUE, pattern = "mHM_restart_001.nc", full.names = TRUE)
  outdir = here("PROC/RAST/calibrated_model_output",m,"AWC")
  dir.create(outdir)
  for (i in 1:n) {
    # i=2
    var = variables[i]
    for (j in 1:length(files)) {
      x = nc_open(files[j]);x
      
      # extraer atributos de las variables
      global.att = ncatt_get(x, 0)
      v.att = ncatt_get(x, var)
      
      # extraer variables
      v = ncvar_get(x, var)
      
      lon = seq(global.att$xllcorner_L1, 
                global.att$xllcorner_L1+((global.att$nrows_L1-1)*global.att$cellsize_L1),
                0.03125) %>% as.array
      
      lat = seq(global.att$yllcorner_L1, 
                global.att$yllcorner_L1+((global.att$ncols_L1-1)*global.att$cellsize_L1),
                0.03125) %>% rev %>% as.array
      if (m == "FAO.calib"){
        soilh = c(30,100) %>% as.array
      }
      else if (m == "SG.calib"){
        soilh = c(5,15,30,60,100) %>% as.array
        
      } else{
        soilh = c(5,15,30,60,100,200) %>% as.array
        
      }
      dim(v);dim(lon);dim(lat);dim(soilh)
      
      # close file
      nc_close(x)
      
      # create and write the netCDF file -- ncdf4 version ---------------------
      
      # define dimensions
      londim <- ncdim_def("lon","degrees_east",as.double(lon)) 
      latdim <- ncdim_def("lat","degrees_north",as.double(lat))
      soildim <- ncdim_def("soilh","bottom_cm", as.double(soilh))
      
      v.att
      # define variables
      fillvalue = v.att[['_FillValue']]
      if (var == "L1_soilMoistFC"){
        v_def <- ncvar_def(name = var,units = 'mm', 
                           dim = list(londim,latdim,soildim),
                           missval = fillvalue,
                           longname = 'Field Capacity', 
                           prec="double",
                           verbose = TRUE)
      }else{
        v_def <- ncvar_def(name = var,units = 'mm', 
                           dim = list(londim,latdim,soildim),
                           missval = fillvalue,
                           longname = 'Permanent Wilting Point', 
                           prec="double",
                           verbose = TRUE)
      }
      
      # create netCDF file and put arrays
      ncfname = paste0(outdir, '/',var,"_",j,'.nc')
      ncout <- nc_create(ncfname,list(v_def),force_v4=TRUE)
      
      # put variables
      ncvar_put(ncout,v_def,v)
      
      # put additional attributes into dimension and data variables
      ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
      ncatt_put(ncout,"lat","axis","Y")
      ncatt_put(ncout,"soilh","axis","H")
      
      # Get a summary of the created file:
      ncout
      
      # close the file, writing data to disk
      nc_close(ncout)
      
    }
  }
}
