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
                "scico") # list of packages to load
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
  conflicted::conflict_prefer("last", "dplyr")
  # LOAD FUNCTIONS ------------------------------------
}
m = "CLSM"
files = list.files(here("proc","rast","calibrated_model_output",m,"aET"), full.names = TRUE, pattern = "tif$")
index = grep("H6", files)
files = files[-index]
length(files)

files.it = files %>% 
  str_split("_") %>%
  sapply(., function(x){x[[5]]}) %>%
  unlist %>% 
  unique

out = here("PROC","RAST","CALIBRATED_MODEL_OUTPUT",m,"aET_0.100")
dir.create(out)
lapply(files.it, function(f){
  index = grep(f, files)
  r1 = rast(files[index[1]])
  r2 = rast(files[index[2]])
  r3 = rast(files[index[3]])
  r4 = rast(files[index[4]])
  r5 = rast(files[index[5]])
  rt = sum(r1,r2,r3,r4,r5)
  writeRaster(rt, filename = here(
    out, 
    str_c("aET_monthly_",f,"_0.100.cm.tif")
  ))
})

files = list.files(here("proc","rast","calibrated_model_output",m,"SM"), full.names = TRUE, pattern = "tif$")
index = grep("H6", files)
files = files[-index]
length(files)

files.it = files %>% 
  str_split("_") %>%
  sapply(., function(x){x[[5]]}) %>%
  unlist %>% 
  unique

out = here("PROC","RAST","CALIBRATED_MODEL_OUTPUT",m,"SWC_0.100")
dir.create(out)
lapply(files.it, function(f){
  index = grep(f, files)
  r1 = rast(files[index[1]])*50
  r2 = rast(files[index[2]])*100
  r3 = rast(files[index[3]])*150
  r4 = rast(files[index[4]])*300
  r5 = rast(files[index[5]])*400
  rt = sum(r1,r2,r3,r4,r5)
  writeRaster(rt, filename = here(
    out, 
    str_c("SWC_monthly_",f,"_0.100.cm.tif")
  ), overwrite=TRUE)
})

