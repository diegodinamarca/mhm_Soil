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
                "viridis",
                "scico"
  ) # list of packages to load
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
  source(here("SCRIPTS","helper_functions_2.R"))
}

# Extract data from fluviometric station ----------------------------------
sites = read_sf(here("DATA","SHP","estaciones_fluvio.geojson")) %>% 
  st_transform(4326) %>% 
  rename(sitename = estacion) %>% 
  select(sitename)
sites

m = "SG"
folders = list.files(here("BASIN_DATA/out/7339001/calibration_runs/",str_c(m,"_calibrated/")), full.names = TRUE)
folders
out = here("PROC","TABLE","calibrated_model_output","estfluv",m)
dir.create(out)

# no procesar outputs que ya fueron procesados anteriormente
reprocess = FALSE
if (!reprocess){
  files.done = lapply(list.files(here(out)), function(f){
    str_split(f,"_") %>% .[[1]] %>% xts::last() %>% str_sub(end = -5)
  }) %>% unlist
  if (is.null(files.done)){
    message("Can't find any processed files, processing all outputs")
  }else{
    index = sapply(files.done, function(f){
      pat = str_c(f, "$")
      grep(pat, folders)
    }) %>% as.vector()
    if (length(index)!= 0){
      folders = folders[-index]
      
    }else{
      message("Index vector lengths = 0. Maybe name of processed files don't match output folder.")
    }
  }
  
}
folders
for (i in 1:length(folders)) {
  
  message("iteration: ",i, "/", length(folders))
  folder = folders[i]
  
  # model iteration
  model.it = xts::last(str_split(folder, "/")[[1]]) %>% 
    str_sub(start = 8)
  
  # extract Qb
  flux = get_flux(folder, m, n = 11, daily = TRUE, process = FALSE)
  df.roff = sites_extract(flux, sites)
  names(df.roff) = c("sitename","time", "roff")
  
  # extract Qf
  flux = get_flux(folder, m, n = 13, daily = TRUE, process = FALSE)
  df.qf =  sites_extract(flux, sites)
  names(df.qf) = c("sitename","time", "Qf")
  
  # extract Qs
  flux = get_flux(folder, m, n = 14, daily = TRUE, process = FALSE)
  df.qs = sites_extract(flux, sites)
  names(df.qs) = c("sitename","time", "Qs")
 
  # extract Qb
  flux = get_flux(folder, m, n = 15, daily = TRUE, process = FALSE)
  df.qb = sites_extract(flux, sites)
  names(df.qb) = c("sitename","time", "Qb")
  
  df = full_join(df.roff, df.qf) %>% 
    full_join(df.qs) %>% 
    full_join(df.qb)
  df$model = model.it
  
  write_csv(df, here(out, str_c("mhm_output_estfluv_",model.it,".csv")))
}


# Model summary -----------------------------------------------------------
in.dir = here("PROC","TABLE","calibrated_model_output","estfluv",m)
in.dir
files = list.files(in.dir, full.names = TRUE)
files
dataset = lapply(files, read_csv) %>% rlist::list.rbind()
dataset

daily_summary = dataset %>%
  pivot_longer(cols = roff:Qb, names_to = "var", values_to = "value") %>% 
  group_by(time, var, sitename) %>% 
  summarise(
    min = min(value),
    q10 = quantile(value, 0.1),
    median = median(value),
    mean = mean(value),
    sd = sd(value),
    q90 = quantile(value, 0.9),
    max = max(value)
  )
daily_summary
beepr::beep(sound = 5)
write_csv(daily_summary, here("PROC","TABLE", str_c("mhm_output_summary_estfluv_",m,".csv")))




