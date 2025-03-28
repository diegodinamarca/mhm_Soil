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

library(soiltexture)

# Directorio de los mapas de suelo input de mHM
dir = here("BASIN_DATA","data","7339001")
model.dir = here(dir, c("morph_cr2soil","morph_soilgrids","morph"))
model.dir
lut.dir = here(dir,str_c("class_definitions_files", c("_cr2soil","_SG","")))
lut.dir
model = c("CLSM","SG","FAO")
out.dir = here("PROC","RAST","TEX.CLASS.MHM",model)
out.dir
sapply(out.dir, dir.create, recursive = TRUE)
m=1
# Tabla con los valores de Arcilla, arena y densidad aparente por ID
for (m in 1:3){
  # m=3
  print(model[m])
  lut.file = here(lut.dir[m],"soil_classdefinition_iFlag_soilDB_1.txt")
  
  lut = read.table(lut.file, skip = 1, header = TRUE) %>% as_tibble
  lut
  names(lut) = c("ID","CLAY","SAND","BD","nSamples")
  
  lut = lut %>% 
    mutate(SILT = abs(100-(SAND+CLAY))) %>% 
    mutate(TOTAL = CLAY+SAND+SILT) %>% 
    as.data.frame()
  lut %>% head
  
  model.files = list.files(
    model.dir[m],
    pattern = "[0-9].asc",
    full.names = TRUE
  )
  
  # i=1
  x = rast(model.files)
  v = x %>% 
    values(dataframe = TRUE, na.rm = FALSE)

  # Transformar a mapas de arcilla
  clay = lut$CLAY
  names(clay) = lut$ID
  r.clay = x
  v.clay = v %>% mutate(across(starts_with("soil"), ~ clay[.x]))
  v.clay
  values(r.clay) = v.clay
  depths.width = c(5,10,15,30,40,100)
  r.clay.sum = round(sum(r.clay*depths.width*(1/200)), digits = 0)
  r.clay.sum.up = round(sum(r.clay[[1:5]]*depths.width[1:5]*(1/100)), digits = 0)
  r.clay.sum.deep = round(sum(r.clay[[6]]*depths.width[6]*(1/100)), digits = 0)
  
  r.clay.sum
  # plot(r.clay.sum)
  
  sand = lut$SAND
  names(sand) = lut$ID
  r.sand = x
  v.sand = v %>% mutate(across(starts_with("soil"), ~ sand[.x]))
  v.sand
  values(r.sand) = v.sand
  depths.width = c(5,10,15,30,40,100)
  r.sand.sum = round(sum(r.sand*depths.width*(1/200)), digits = 0)
  r.sand.sum.up = round(sum(r.sand[[1:5]]*depths.width[1:5]*(1/100)), digits = 0)
  r.sand.sum.deep = round(sum(r.sand[[6]]*depths.width[6]*(1/100)), digits = 0)
  
  r.sand.sum
  # plot(r.sand.sum)
  
  silt = lut$SILT
  names(silt) = lut$SILT
  # silt[v$soil_class_horizon_01]
  r.silt = x
  v.silt = v %>% mutate(across(starts_with("soil"), ~ silt[.x]))
  v.silt
  values(r.silt) = v.silt
  depths.width = c(5,10,15,30,40,100)
  r.silt.sum = round(sum(r.silt*depths.width*(1/200)), digits = 0)
  r.silt.sum.up = round(sum(r.silt[[1:5]]*depths.width[1:5]*(1/100)), digits = 0)
  r.silt.sum.deep = round(sum(r.silt[[6]]*depths.width[6]*(1/100)), digits = 0)

  r.silt.sum
  # plot(r.silt.sum)
  

  classify_texture <- function(r) {
    names(r) = c("CLAY","SAND","SILT")
    plot(r)
    writeRaster(r, here(out.dir[m], "horizons_weighted_mean.tif"), overwrite = TRUE)
    
    tex.data = values(r, dataframe = TRUE, na.rm = TRUE)
    tex.norm = TT.normalise.sum(tri.data = tex.data[c('CLAY', 'SILT',  'SAND')])  %>% 
      rename(SAND_n = SAND, CLAY_n = CLAY, SILT_n = SILT)
    
    tex.class = TT.points.in.classes(tri.data  = tex.norm,
                                     css.names = c('CLAY_n', 'SILT_n', 'SAND_n'),
                                     class.sys = "USDA.TT", 
                                     PiC.type  = "t",
                                     collapse  = ', ')
    tex.class = sapply(str_split(tex.class,", "), function(vector){vector[1]})
    
    tex.data = bind_cols(tex.data, class = tex.class) %>% 
      mutate(class.id = as.integer(factor(class)))
    
    
    r.class = raster(r)
    v.na = values(r[[1]], dataframe = FALSE, na.rm = FALSE)
    index = which(!is.na(v.na))
    r.class[index] = tex.data$class.id
    return(list(rclass = r.class, texdata = tex.data))
  }
  
  r = c(r.clay.sum, r.sand.sum, r.silt.sum)
  r.up = c(r.clay.sum.up, r.sand.sum.up, r.silt.sum.up)
  r.deep = c(r.clay.sum.deep, r.sand.sum.deep, r.silt.sum.deep)
  
  res = classify_texture(r)
  writeRaster(res$rclass, here(out.dir[m], "horizons_weighted_mean_tex.tif"), overwrite = TRUE)
  lut.out = res$texdata %>% 
    group_by(class, class.id) %>% 
    slice(1) %>% 
    select(class, class.id)
  write_csv(lut.out, here(out.dir[m],"LUT.csv"))
  
  res = classify_texture(r.up)
  writeRaster(res$rclass, here(out.dir[m], "horizons_weighted_mean_tex_up.tif"), overwrite = TRUE)
  lut.out = res$texdata %>% 
    group_by(class, class.id) %>% 
    slice(1) %>% 
    select(class, class.id)
  write_csv(lut.out, here(out.dir[m],"LUT_up.csv"))
  
  res = classify_texture(r.deep)
  writeRaster(res$rclass, here(out.dir[m], "horizons_weighted_mean_tex_deep.tif"), overwrite = TRUE)
  
  lut.out = res$texdata %>% 
    group_by(class, class.id) %>% 
    slice(1) %>% 
    select(class, class.id)
  write_csv(lut.out, here(out.dir[m],"LUT_deep.csv"))

}
  





