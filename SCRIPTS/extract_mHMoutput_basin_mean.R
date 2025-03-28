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

get_flux = function(folder, model, n = 4, daily, process = TRUE) {
  l = extract_model_q(folder, model)
  y.start = l$time_start
  y.end = l$time_end
  if (daily){
    yearmonth = seq(y.start, y.end, by="days")
  }else{
    yearmonth = seq(y.start, y.end, by="months")
  }
  # print("yearmonth")
  cat(paste0("processin outputflux_",n,"\n"))
  f = do.call(paste0("get_outFstate_",n), list(folder, yearmonth))
  return(f)
  # f = get_outFstate_(model, idbasin, yearmonth, outdir)
}

basin_mean <- function(flux) {
  as_tibble(values(flux)) %>% 
    pivot_longer(cols = 1:ncol(.), names_to = "time", values_to = "var") %>% 
    mutate(time = ymd(time)) %>%
    drop_na() %>%
    group_by(time) %>%
    summarise(var = mean(var))
}

basin_mean_bylayer <- function(flux) {
  data_day = lapply(1:6, function(i){
    if (!is.null(flux[[i]])){
      x = as_tibble(values(flux[[i]])) %>% 
        pivot_longer(cols = 1:nlyr(flux[[i]]), names_to = "time",values_to = "var") %>% 
        drop_na() %>% 
        mutate(time = ymd(time)) %>% 
        group_by(time) %>%
        summarise(var = mean(var))
    }
  }) %>% bind_cols()
  data_day = data_day %>% 
    dplyr::select(time = 1, starts_with("var"))
}

# Modelo a procesar
m = "SG"

# Carpeta con los archivos de salida de mHM
folders = list.files(here("BASIN_DATA/out/7339001/calibration_runs/",str_c(m,"_calibrated/")), full.names = TRUE)
folders

# Carpeta de salida
out = here("PROC","TABLE","calibrated_model_output","basin_mean",m)
dir.create(out)

# Cuenca ------------------------------------------------------------------
cuenca = read_sf("DATA/SHP/Cauquenes.shp") %>% st_transform(4326)


# Procesa forzantes del modelo. Solo hay que hacerlo una vez a menos que cambien las forzantes
process.forcings=FALSE


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
# Extract aET and SWC by layer and calculate mean of all pixels (basin mean)
for (i in 1:length(folders)) {
  message("iteration: ",i, "/", length(folders))
  folder = folders[i]
  
  # model iteration
  model.it = xts::last(str_split(folder, "/")[[1]]) %>% 
    str_sub(start = 8)
  
  # extract aET
  flux = get_flux(folder, m, n = 19, daily = TRUE, process = FALSE)
  df.et = basin_mean_bylayer(flux)
  names(df.et)[2:7] = paste0("aET_L", 1:6)
  
  # extract SWC
  flux = get_flux(folder, m, n = 3, daily = TRUE, process = FALSE)
  df.swc = basin_mean_bylayer(flux)
  names(df.swc)[2:7] = paste0("swc_L", 1:6)
  
  df = full_join(df.et, df.swc)
  df$model = model.it
  
  write_csv(df, here(out, str_c("mhm_output_basin_mean_",model.it,".csv")))
}

# Forzantes ---------------------------------------------------------------
if (process.forcings){
  #   Precipitacion -----------------------------------------------------------
  pp = rast(here("BASIN_DATA","DATA","7339001","meteo_cr2met_v2.5_1960_2023","pre.NC"))
  names(pp) = time(pp)
  
  pp_day = extract(pp, cuenca, fun = "mean") %>% 
    as_tibble() %>% 
    pivot_longer(cols = 2:ncol(.), names_to = "time",values_to = "PP") %>% 
    mutate(time = ymd(time)) %>% 
    select(-ID)
  pp_day
  data_day = pp_day
  
  #   PET ---------------------------------------------------------------------
  dates <- seq(ymd("1960-01-01"),ymd("2021-12-31"), by = "days");length(dates)
  pet = rast("D:/mHM/BASIN_DATA/data/7339001/meteo_cr2met_v2.5_1960_2023/pet.nc")
  names(pet) <- time(pet)
  pet_day = extract(pet, cuenca, fun = "mean") %>% 
    as_tibble() %>% 
    pivot_longer(cols = 2:ncol(.), names_to = "time",values_to = "PET_cr2") %>% 
    mutate(time = ymd(time)) %>% 
    select(-ID)
  pet_day
  
  data_day = full_join(data_day, pet_day, by = c("time"))
  data_day
  write_csv(data_day, here("PROC","TABLE","mhm_forcings_basin_mean_daily.csv"))
  
}


# model statistics --------------------------------------------------------
m = "SG"
in.dir = here("PROC","TABLE","calibrated_model_output","basin_mean",m)

files = list.files(in.dir, full.names = TRUE)
files
dataset = lapply(files, read_csv) %>% rlist::list.rbind()
dataset

daily_summary = dataset %>%
  pivot_longer(cols = aET_L1:swc_L6, names_to = "var", values_to = "value") %>% 
  group_by(time, var) %>% 
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

write_csv(daily_summary, here("PROC","TABLE", str_c("mhm_output_summary_basin_mean_",m,".csv")))

daily_summary_2 = dataset %>%
  mutate(aET_L1 = aET_L1+aET_L2+aET_L3+aET_L4+aET_L5,
         aET_L2 = aET_L6,
         swc_L1 = swc_L1+swc_L2+swc_L3+swc_L4+swc_L5,
         swc_L2 = swc_L6) %>%
  
  select(time, model, aET_L1, aET_L2, swc_L1, swc_L2) %>% 
  pivot_longer(cols = aET_L1:swc_L2, names_to = "var", values_to = "value") %>% 
  group_by(time, var) %>% 
  summarise(
    min = min(value),
    q10 = quantile(value, 0.1),
    median = median(value),
    mean = mean(value),
    sd = sd(value),
    q90 = quantile(value, 0.9),
    max = max(value)
  )
write_csv(daily_summary_2, here("PROC","TABLE", str_c("mhm_output_summary_basin_mean_",m,"_2layers.csv")))

df.iterations = dataset %>%
  mutate(aET_L1 = aET_L1+aET_L2+aET_L3+aET_L4+aET_L5,
         aET_L2 = aET_L6,
         swc_L1 = swc_L1+swc_L2+swc_L3+swc_L4+swc_L5,
         swc_L2 = swc_L6) %>%
  
  select(time, model, aET_L1, aET_L2, swc_L1, swc_L2) %>% 
  # pivot_longer(cols = aET_L1:swc_L2, names_to = "var", values_to = "value") %>% 
  mutate(time = floor_date(time, unit = "months")) %>% 
  group_by(time, model) %>% 
  summarise(across(starts_with("aET"), sum),
            across(starts_with("swc"), mean)) %>% 
  mutate(model.it = model,
         model = m)
df.iterations
write_csv(df.iterations, here("PROC","TABLE", str_c("mhm_output_basin_mean_monthly_",m,"_2layers.csv")))


# aET plots ---------------------------------------------------------------
daily_summary = read_csv(here("PROC","TABLE", str_c("mhm_output_summary_basin_mean_",m,".csv")))

monthly_sum = daily_summary %>% 
  filter(str_detect(var, "aET")) %>% 
  mutate(time_month = floor_date(time, unit = "month")) %>% 
  group_by(time_month, var) %>% 
  summarise(across(where(is.numeric), sum))

monthly_sum %>% 
  ggplot(aes(x = time_month))+
  facet_wrap(.~var, scales = "free_y")+
  geom_line(aes(y = mean))+
  geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = 0.5)+
  scale_x_date(limits = c(ymd("2020-01-01"),ymd("2023-12-01")))


monthly_ltm = monthly_sum %>% 
  mutate(m = month(time_month)) %>% 
  group_by(m, var) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  mutate(m = factor(m, levels = c(5:12, 1:4)))

ggplot(monthly_ltm, aes(x = as.numeric(m)))+
  facet_wrap(.~var, scales = "free_y")+
  geom_line(aes(y = mean))+
  geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = 0.5)+
  scale_x_continuous(breaks = 1:12, labels = month.abb[c(5:12, 1:4)])+
  labs(x= "",y = "aET [mm]", title = "LTM monthly aET of CLSM model across 90 calibration runs",
       caption = "The ribbon represents the [mean-sd, mean+sd] range")


yearly_sum = monthly_sum %>% 
  mutate(time_year = floor_date(time_month, unit = "year")) %>% 
  group_by(time_year, var) %>% 
  summarise(across(where(is.numeric), sum))

yearly_sum %>% 
  ggplot(aes(x = time_year))+
  facet_wrap(.~var, scales = "free_y")+
  geom_line(aes(y = mean))+
  geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = 0.5)+
  scale_x_date(limits = c(ymd("1965-01-01"),ymd("2023-12-01")))+
  labs(x= "",y = "aET [mm]", title = "LTM monthly aET of CLSM model across 90 calibration runs",
       caption = "The ribbon represents the [mean-sd, mean+sd] range")


yearly_sum %>%
  group_by(time_year) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  ggplot(aes(x = time_year))+
  # facet_wrap(.~var, scales = "free_y")+
  geom_line(aes(y = mean))+
  geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = 0.5)+
  scale_x_date(limits = c(ymd("1965-01-01"),ymd("2023-12-01")))+
  labs(x= "",y = "aET [mm]", title = "Annual aET of CLSM model across 90 calibration runs",
       caption = "The ribbon represents the [mean-sd, mean+sd] range")


# swc plots ---------------------------------------------------------------

monthly_sum = daily_summary %>% 
  filter(str_detect(var, "swc")) %>% 
  mutate(time_month = floor_date(time, unit = "month")) %>% 
  group_by(time_month, var) %>% 
  summarise(across(where(is.numeric), mean))

monthly_sum %>% 
  ggplot(aes(x = time_month))+
  facet_wrap(.~var, scales = "free_y")+
  geom_line(aes(y = mean))+
  geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = 0.5)+
  scale_x_date(limits = c(ymd("2020-01-01"),ymd("2023-12-01")))

monthly_ltm = monthly_sum %>% 
  mutate(m = month(time_month)) %>% 
  group_by(m, var) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  mutate(m = factor(m, levels = c(5:12, 1:4)))

ggplot(monthly_ltm, aes(x = as.numeric(m)))+
  facet_wrap(.~var, scales = "free_y")+
  geom_line(aes(y = mean))+
  geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = 0.5)+
  scale_x_continuous(breaks = 1:12, labels = month.abb[c(5:12, 1:4)])+
  labs(x= "",y = "SWC [mm]", title = "LTM monthly SWC of CLSM model across 90 calibration runs",
       caption = "The ribbon represents the [mean-sd, mean+sd] range")

yearly_sum = monthly_sum %>% 
  mutate(time_year = floor_date(time_month, unit = "year")) %>% 
  group_by(time_year, var) %>% 
  summarise(across(where(is.numeric), mean))

yearly_sum %>% 
  ggplot(aes(x = time_year))+
  facet_wrap(.~var, scales = "free_y")+
  geom_line(aes(y = mean))+
  geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = 0.5)+
  scale_x_date(limits = c(ymd("1965-01-01"),ymd("2023-12-01")))+
  labs(x= "",y = "SWC [mm]", title = "Annual mean SWC of CLSM model across 90 calibration runs",
       caption = "The ribbon represents the [mean-sd, mean+sd] range")


yearly_sum %>%
  group_by(time_year) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  ggplot(aes(x = time_year))+
  # facet_wrap(.~var, scales = "free_y")+
  geom_line(aes(y = mean))+
  geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = 0.5)+
  scale_x_date(limits = c(ymd("1965-01-01"),ymd("2023-12-01")))+
  labs(x= "",y = "swc [mm]", title = "Annual mean SWC of CLSM model across 90 calibration",
       caption = "The ribbon represents the [mean-sd, mean+sd] range")
