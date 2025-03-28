library(terra)
library(tidyverse)
library(sf)
library(lubridate)
library(hydroTSM)
library(sf)
library(hydroGOF)
library(scico)
library(rasterVis)

extract_sim_streamflow_iterations <- function(folder) {
  files = list.files(folder, 
                     full.names = TRUE, recursive = TRUE, pattern = "daily_discharge.out")
  files
  for (i in 1:length(files)) {
    f = files[i]
    x = read_table(f)
    x$date = as.Date(paste(x$Day, x$Mon, x$Year, sep = " "), format = "%d %m %Y")
    x = tibble(date = x$date, qobs = x[[5]], qsim = x[[6]])
    # setting NA valuss
    misVal = -9999.0
    
    if (i == 1){
      qday = x %>% mutate(qobs = ifelse(qobs == -9999, NA, qobs),
                          qsim = ifelse(qsim == -9999, NA, qsim))
      df = qday
    }else{
      qday = x %>% mutate(qobs = ifelse(qobs == -9999, NA, qobs),
                          qsim = ifelse(qsim == -9999, NA, qsim)) %>% 
        select(date, qsim)
      df = full_join(df, qday, by = "date")
      
    }
  }
  names(df)[3:nrow(df)] = str_c("qsim_", 1:(nrow(df)-2))
  df
}

# Function to calculate error metrics
calculate_error_metrics <- function(obs, sim) {
  # browser()
  N <- length(obs)
  bias <- mean(sim - obs)
  pbias <- hydroGOF::pbias(sim, obs)
  rmse <- rmse(sim, obs)
  correlation <- cor(sim, obs)
  ubrmse <- sqrt(mean((sim - mean(sim) - (obs - mean(obs)))^2))
  nse <- NSE(sim, obs)
  kge_components <- KGE(sim, obs, method="2009", out.type = "full")
  kge <- kge_components$KGE.value
  alpha <- kge_components$KGE.elements[["Alpha"]]
  beta <- kge_components$KGE.elements[["Beta"]]
  r <- kge_components$KGE.elements[["r"]]
  hbe <- (sum(sim)-sum(obs))/sum(obs)
  return(tibble(
    N = N,
    bias = bias,
    rmse = rmse,
    correlation = correlation,
    ubrmse = ubrmse,
    nse = nse,
    kge = kge,
    alpha = alpha,
    beta = beta,
    r = r,
    hbe = hbe,
    pbias = pbias
  ))
}

extract_model_q = function(folder, model, plot.dwdata = FALSE){
  file = list.files(folder, pattern = "daily_discharge.out", full.names = TRUE)
  x = read_table(file)
  
  qday = x %>% 
    mutate(date = as.Date(paste(Day, Mon, Year, sep = " "), format = "%d %m %Y")) %>% 
    select(date, qobs = starts_with("Qobs"), qsim = starts_with("Qsim")) %>% 
    mutate(qobs = ifelse(qobs == -9999, NA, qobs),
           qsim = ifelse(qsim == -9999, NA, qsim)) %>% 
    as.data.frame
  
  # diario a anual
  qyear = daily2annual(qday, FUN = mean, na.rm = TRUE)
  qmonth = daily2monthly(qday, FUN = mean, na.rm = TRUE)
  qmonth
  x2 = qday %>% mutate(date = date %m-% months(3))
  qwyear = daily2annual(x2, FUN = mean, na.rm = TRUE)
  
  if (plot.dwdata){
    ndays = x2 %>% tibble %>% 
      mutate(date = floor_date(date, unit = "year")) %>% 
      group_by(date) %>%
      drop_na() %>% 
      summarise(n = n()) %>% 
      mutate(ly = leap_year(date),
             days = if_else(ly, 366, 365),
             prop = 100*n/days)
    write_csv(ndays, str_c("PROC/days_with_data_",idbasin,"_wateryear.csv"))
    ndays %>% 
      # filter(prop > 80) %>%
      ggplot()+
      geom_line(aes(x = date, y = prop))+
      geom_point(aes(x = date, y = prop))+
      labs(title = str_c("Proportion of days with data, basin ID: ",idbasin))+
      scale_x_date(date_breaks = "2 years", date_labels = "%Y")+
      scale_y_continuous(n.breaks = 10)
    # theme(axis.text.x = element_text(angle = 90, vjust= 0.5))
    ggsave(str_c("RESULTADOS/FIGS/days_with_data_",idbasin,"_wateryear.png"))
  }
  dates.y = qday$date %>% year() %>% unique
  dates.wy = qday$date %>% year() %>% unique
  y.start = dates.y[1];y.start
  y.end = dates.y[length(dates.y)];y.end
  dates.m = qday$date %>% floor_date(unit = "month") %>% unique
  list(day = tibble(model= model, qday), 
       month = qmonth,
       year = qyear,
       wateryear = qwyear,
       time_start = ymd(qday$date %>% dplyr::first()),
       time_end = ymd(qday$date %>% dplyr::last()))
}

# Actual Evapotransporation from the soil layers
get_outFstate_19 = function(folder, yearmonth = NULL){
  flux.nam = "aET"
  layers = c("L01","L02","L03","L04","L05","L06")
  
  files = list.files(folder,full.names = TRUE, pattern = ".nc$");files
  # read fluxes
  index = grep("mHM_Fluxes_States", files)
  fluxes = rast(files[index])
  fluxes
  # names(fluxes) %>% as.data.frame %>% write.csv("nombre flujos.csv")
  
  rl = lapply(layers, function(layer_i){
    index = grep(paste0(flux.nam,"_",layer_i), names(fluxes))
    r = NULL
    if (length(index) != 0){
      r = fluxes[[index]]
      tryCatch({names(r) = yearmonth}, 
               error = function(e){
                 message("Fecha de inicio y/o fin no coinciden con el flujo de mHM")}
      )
    }else{
      cat("No se encuentra ", paste0(flux.nam,"_",layer_i), " entre los outputs de la cuenca ", 
          idbasin, "para el modelo ", model,"\n")
    }
    r
  })
  names(rl) = paste0(flux.nam,"_",layers)
  return(rl)
  
  
}
# baseflow generated per cell        (L1_baseflow)     [mm/T] 
get_outFstate_15 = function(folder, yearmonth = NULL){
  
  flux.nam = "QB"
  files = list.files(folder,full.names = TRUE, pattern = ".nc$");files
  
  # read fluxes
  index = grep("mHM_Fluxes_States", files)
  fluxes = rast(files[index])
  fluxes
  # names(fluxes) %>% as.data.frame %>% write.csv("nombre flujos.csv")
  
  index = grep(flux.nam, names(fluxes))
  r = NULL
  if (length(index) != 0){
    r = fluxes[[index]]
    tryCatch({names(r) = yearmonth}, 
             error = function(e){
               message("Fecha de inicio y/o fin no coinciden con el flujo de mHM")}
    )
  }else{
    cat("No se encuentra ", flux.nam, " entre los outputs de la cuenca ",
        idbasin, "para el modelo ", model,"\n")
  }
  
  return(r)
}
# slow interflow generated per cell  (L1_slowRunoff)   [mm/T]
get_outFstate_14 = function(folder, yearmonth = NULL){
  
  flux.nam = "QIs_"
  files = list.files(folder,full.names = TRUE, pattern = ".nc$");files
  
  # read fluxes
  index = grep("mHM_Fluxes_States", files)
  fluxes = rast(files[index])
  fluxes
  # names(fluxes) %>% as.data.frame %>% write.csv("nombre flujos.csv")
  
  index = grep(flux.nam, names(fluxes))
  r = NULL
  if (length(index) != 0){
    r = fluxes[[index]]
    tryCatch({names(r) = yearmonth}, 
             error = function(e){
               message("Fecha de inicio y/o fin no coinciden con el flujo de mHM")}
    )
  }else{
    cat("No se encuentra ", flux.nam, " entre los outputs de la cuenca ",
        idbasin, "para el modelo ", model,"\n")
  }
  
  return(r)
}
# fast interflow generated per cell  (L1_fastRunoff)   [mm/T]
get_outFstate_13 = function(folder, yearmonth = NULL){
  
  flux.nam = "QIf_"
  files = list.files(folder,full.names = TRUE, pattern = ".nc$");files
  
  # read fluxes
  index = grep("mHM_Fluxes_States", files)
  fluxes = rast(files[index])
  fluxes
  # names(fluxes) %>% as.data.frame %>% write.csv("nombre flujos.csv")
  
  index = grep(flux.nam, names(fluxes))
  r = NULL
  if (length(index) != 0){
    r = fluxes[[index]]
    tryCatch({names(r) = yearmonth}, 
             error = function(e){
               message("Fecha de inicio y/o fin no coinciden con el flujo de mHM")}
    )
  }else{
    cat("No se encuentra ", flux.nam, " entre los outputs de la cuenca ",
        idbasin, "para el modelo ", model,"\n")
  }
  
  return(r)
}
# Total discharge generated per cell (L1_total_runoff) [mm/T]
get_outFstate_11 = function(folder, yearmonth = NULL){
  flux.nam = "Q_"
  files = list.files(folder,full.names = TRUE, pattern = ".nc$");files
  
  # read fluxes
  index = grep("mHM_Fluxes_States", files)
  fluxes = rast(files[index])
  fluxes
  # names(fluxes) %>% as.data.frame %>% write.csv("nombre flujos.csv")
  
  index = grep(flux.nam, names(fluxes))
  r = NULL
  if (length(index) != 0){
    r = fluxes[[index]]
    tryCatch({names(r) = yearmonth}, 
             error = function(e){
               message("Fecha de inicio y/o fin no coinciden con el flujo de mHM")}
    )
    
  }else{
    cat("No se encuentra ", flux.nam, " entre los outputs de la cuenca ",
        idbasin, "para el modelo ", model,"\n")
  }
  
  return(r)
}

# ! soil water content in the single layers     (L1_soilMoist)         -- case  3
get_outFstate_3 = function(folder, yearmonth = NULL){
  flux.nam = "SWC"
  layers = c("L01","L02","L03","L04","L05","L06")
  files = list.files(folder,full.names = TRUE, pattern = ".nc$");files
  
  # read fluxes
  index = grep("mHM_Fluxes_States", files)
  fluxes = rast(files[index])
  fluxes
  # names(fluxes) %>% as.data.frame %>% write.csv("nombre flujos.csv")
  
  rl = lapply(layers, function(layer_i){
    index = grep(paste0(flux.nam,"_",layer_i), names(fluxes))
    r = NULL
    if (length(index) != 0){
      r = fluxes[[index]]
      tryCatch({names(r) = yearmonth}, 
               error = function(e){
                 message("Fecha de inicio y/o fin no coinciden con el flujo de mHM")}
      )
    }else{
      cat("No se encuentra ", flux.nam, " entre los outputs de la cuenca ", 
          idbasin, "para el modelo ", model,"\n")
    }
    r
  })
  names(rl) = paste0(flux.nam,"_",layers)
  return(rl)
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

sites_extract <- function(flux, sites) {
  as_tibble(extract(flux, vect(sites))) %>% 
    mutate(sitename = sites$sitename) %>%
    select(ID, sitename, everything()) %>% 
    pivot_longer(cols = 3:ncol(.), names_to = "time",values_to = "var") %>% 
    mutate(time = ymd(time)) %>% 
    select(-ID)
}

sites_extract_bylayer <- function(flux, sites) {
  extr.list = lapply(flux, function(f){
    x = NULL
    if (!is.null(f)){
      x = terra::extract(f, vect(sites))
    }
    x
  })
  data_day = lapply(1:6, function(i){
    if (!is.null(extr.list[[i]])){
      as_tibble(extr.list[[i]]) %>% 
        mutate(sitename = sites$sitename) %>% 
        pivot_longer(cols = 2:ncol(extr.list[[i]]), names_to = "time",
                     values_to = "var") %>% 
        mutate(time = ymd(time)) %>% 
        select(-ID)
    }
  }) %>% bind_cols()
  data_day = data_day %>% 
    select(sitename = 1, time = 2, starts_with("var"))
  data_day
}

  # Funcion para pasar los flujos de diarios a mensuales, por defecto calcular la
# suma mensual
daily_to_monthly = function(r, dates.seq, fun = "sum"){
  ym = NULL
  if (!is.null(r)){
    names(r) = dates.seq
    ym=lapply(year(dates.seq) %>% unique, function(y){
      # ym=lapply(1960:1961, function(y){ 
      # y = 1960
      print(paste0("year ", y))
      index_y = grep(y, year(dates.seq))
      ry = r[[index_y]]
      # cat("month ")
      rym = lapply(month(names(ry)) %>% unique, function(m){
        # cat(m, " ")
        index_m = grep(paste0("^",m,"$"), month(names(ry)))
        if (fun == "sum"){
          rm = ry[[index_m]] %>% sum
        }else if (fun == "mean"){
          rm = ry[[index_m]] %>% mean
        }else if (fun == "sd"){
          rm = ry[[index_m]] %>% stdev
        }else{
          cat("Especificar argumento fun = \"mean\" o \"sum\" o \"sd\" \n")
        }
      }) %>% rast
    }) %>% rast
  }
  
  return(ym)
}

# Funcion para calcular el promedio y acumulado anual de una serie de raster mensuales
# Calcula la suma y el promedio anual de una serie de tiempo de rasters
# years.seq es un vector con los años a procesar
mean_sum_anual = function(r, years.seq){

  r.sum = list()
  r.mean = list()
  for (i in 1:length(years.seq)) {
    y = years.seq[i]
    index = grep(paste0("^",y,"$"), year(names(r)))
    sum.anual = r[[index]] %>% sum
    r.sum[[i]] = sum.anual
    
    mean.anual = r[[index]] %>% mean
    r.mean[[i]] = mean.anual
    
  }
  r.sum = rast(r.sum)
  r.mean = rast(r.mean)
  
  names(r.sum) = years.seq
  names(r.mean) = years.seq
  return(list(anual_sum = r.sum, anual_mean = r.mean))
}

aggregate_sitename = function(sm_sondas){
  sm_sondas %>% 
    mutate(sitename = if_else(sitename == "BN_MOLCO1","BN_MOLCO",
                              if_else(sitename == "BN_MOLCO2", "BN_MOLCO",
                                      if_else(sitename == "BN_MOLCO3", "BN_MOLCO",
                                              if_else(sitename == "PI_MOLCO1", "PI_MOLCO",
                                                      if_else(sitename == "PI_MOLCO1", "PI_MOLCO",
                                                              if_else(sitename == "PI_MOLCO2","PI_MOLCO",
                                                                      if_else(sitename == "MAT_SASB1","MAT_SASB",
                                                                              if_else(sitename == "MAT_SASB2", "MAT_SASB", 
                                                                                      if_else(sitename == "PI_SASB1","PI_SASB",
                                                                                              if_else(sitename == "PI_SASB2","PI_SASB",
                                                                                                      sitename)))))))))))
}

match_standard_horizon=function(x){
  if (x <= 5){
    return(1)
  }else if(x <= 15){
    return(2)
  }else if(x <= 30){
    return(3)
  }else if(x <= 60){
    return(4)
  }else if(x <= 100){
    return(5)
  }else if(x <= 200){
    return(6)
  }else if(x > 200){
    return(6)
  }
}

to_monthly_mean = function(x, dates){
  # obtener solo el mes de cada fecha
  meses = month(as_date(dates))
  
  # raster vacio para guardar resultados
  r = rast()
  for (m in 1:12) {
    index = grep(pattern = str_c("^",m,"$"), meses)
    xm = mean(x[[index]])
    r = c(r, xm)
  }
  names(r) = c("ENE","FEB","MAR","ABR","MAY","JUN","JUL","AGO","SEP","OCT","NOV","DIC")
  return(r)
}

# Funcion para calcular imagenes anuales a partir de imagenes mensuales
# x : imagenes mensuales
# dates : vector con las fechas de cada imagen en formato YYYY-MM-DD
# fun : funcion a aplicar a las imagenes diarias "mean" para media, "sum" para suma
to_wateryear = function(x, dates, fun = "mean"){
  # aproxima al primer dia del mes de cada fecha en dates
  y = floor_date(as_date(dates) %m-% months(3), unit = "year")
  
  # lista de meses unicos
  lista.y = y %>% unique
  
  # numero de meses
  n = length(lista.y)
  
  # dummy raster para guardar los valores mensuales
  x.anual = rast()
  
  # iterar i de 1 a n => i = {1,2,3,...n}
  for (i in 1:n) {
    # obtener la posicion de las imagenes que coinciden con el mes_i
    posicion = grep(lista.y[i], y)
    # calcular la media de todos las imagenes del mes_i
    if (fun == "mean"){
      r = mean(x[[posicion]], na.rm = TRUE)
    }else if(fun == "sum"){
      r = sum(x[[posicion]], na.rm = TRUE)
    }else{
      message("No se reconoce la funcion utilizada, escriba: fun = sum o mean como argumento")
    }
    # agregar el resultado al dummy raster
    x.anual = c(x.anual, r)
  }
  # asignar la fecha a cada imagen
  names(x.anual) = lista.y
  # retornar rasters de valores mensuales
  return(x.anual)
}
to_yearly = function(x, dates, fun = "mean"){
  # aproxima al primer dia del mes de cada fecha en dates
  y = floor_date(as_date(dates), unit = "year")
  
  # lista de meses unicos
  lista.y = y %>% unique
  
  # numero de meses
  n = length(lista.y)
  
  # dummy raster para guardar los valores mensuales
  x.anual = rast()
  
  # iterar i de 1 a n => i = {1,2,3,...n}
  for (i in 1:n) {
    # obtener la posicion de las imagenes que coinciden con el mes_i
    posicion = grep(lista.y[i], y)
    # calcular la media de todos las imagenes del mes_i
    if (fun == "mean"){
      r = mean(x[[posicion]], na.rm = TRUE)
    }else if(fun == "sum"){
      r = sum(x[[posicion]], na.rm = TRUE)
    }else{
      message("No se reconoce la funcion utilizada, escriba: fun = sum o mean como argumento")
    }
    # agregar el resultado al dummy raster
    x.anual = c(x.anual, r)
  }
  # asignar la fecha a cada imagen
  names(x.anual) = lista.y
  # retornar rasters de valores mensuales
  return(x.anual)
}

reduce_mean <- function(flux, varname, model) {
  n = str_c(varname,"_",model)
  as_tibble(values(flux)) %>% 
    pivot_longer(cols = 1:ncol(.), names_to = "time",values_to = varname) %>% 
    mutate(time = ymd(time)) %>%
    drop_na() %>%
    group_by(time) %>%
    summarise("{varname}_{model}" := mean(get(varname)))
}

mmean_to_seasons = function(r){
  summ = mean(r[[1:3]])
  fall = mean(r[[4:6]])
  wint = mean(r[[7:9]])
  sprin = mean(r[[10:12]])
  res = c(fall, wint, sprin, summ)
  names(res) = c("Fall","Winter","Spring","Summer")
  return(res)
}

mean_from_period <- function(x, start, end) {
  start_y = grep(start,names(x))
  end_y = grep(end, names(x))
  mhm_meany = mean(x[[start_y:end_y]])
}

annotate_metrics <- function(sim, obs) {
  m <- hydroGOF::gof(sim, obs, norm = "maxmin")
  msg <- paste0("R2 = ", m["R2", ], "\nRMSE = ", m["RMSE", ], "\nNRMSE = ", m["NRMSE %", ], "%", "\nPBIAS = ", m["PBIAS %", ], "%")
  as.character(msg)
}
get_season <- function(x) {
  if (x %in% c(1,2,3)){
    s = "summer"
  }else if(x %in% c(4,5,6)){
    s = "fall"
  }else if(x %in% c(7,8,9)){
    s = "winter"
  }else{
    s = "spring"
  }
}
get_season_v = Vectorize(get_season)

daily_climatology = function(x){
  rmean = rast()
  rsd = rast()
  rq.025 = rast()
  rq.975 = rast()
  for (j in yday) {
    message("processing day:", j)
    index = grep(str_c("^",j,"$"), yday(dates))
    xx = x[[index]]
    m = mean(xx)
    rmean = c(rmean, m)
    s = stdev(xx)
    rsd = c(rsd, s)
    q.025 = quantile(xx, 0.025)
    rq.025 = c(rq.025, q.025)
    q.975 = quantile(xx, 0.975)
    rq.975 = c(rq.975, q.975)
    
  }
  names(rmean) = yday
  names(rsd) = yday
  names(rq.025) = yday
  names(rq.975) = yday
  return(list(mean = rmean, sd = rsd,
              q.025 = rq.025, q.975 = rq.975))
}

calc_KGE = function(sim, obs, qmean = NA, daily = TRUE){
  # mKGE = mean( hydroGOF::KGE(sim, obs, method = "2009"), hydroGOF::KGE(1/(sim+qmean), 1/(obs+qmean), method = "2009") )
  # browser()
  # qmean = as.numeric(qmean)
  if (daily){
    message("Calculating daily metrics")
    KGE_1Q = hydroGOF::KGE(1/(sim+qmean), 1/(obs+qmean), method = "2009")
  } else {
    message("Calculating monthly metrics")
    KGE_1Q = hydroGOF::KGE(1/(sim), 1/(obs), method = "2009")
  }
  KGE_Q = hydroGOF::KGE(sim, obs, method = "2009")
  mKGE = mean(c(KGE_Q, KGE_1Q))
  logKGE = hydroGOF::KGE(log(sim+qmean), log(obs+qmean), method = "2009")
  KGE = hydroGOF::KGE(sim, obs, method = "2009")
  NSE = hydroGOF::NSE(sim, obs)
  logNSE = hydroGOF::NSE(log(sim+qmean), log(obs+qmean))
  data.frame(KGE, logKGE, KGE_Q, KGE_1Q, mKGE, NSE, logNSE)
}

clean_names <- function(x) {
  oldnames = x %>% select(ends_with(".x")) %>% names()
  if (length(oldnames) != 0){
    
    newnames = oldnames %>% 
      str_sub(end = -3)
    names(oldnames) = newnames
    
    x = x %>% select(
      -ends_with(".y")
    ) %>% 
      rename(!!!oldnames)
    
    oldnames = names(x)[3:ncol(x)]
    oldnames
    newnames = oldnames %>% 
      str_replace("_","-") %>% 
      str_replace("_","-") 
    
    names(oldnames) = newnames
    
    x = x %>% select(
      -ends_with(".y")
    ) %>% 
      rename(!!!oldnames)
    message("columns ending with .x, removed")
  }else{
    message("names already clean")
  }
  return(x)
}

hour_to_day <- function(x, fun = "mean", mode = "8a8") {
  ts_to_mean_day8a8 = function(x){
    x %>% 
      mutate(
        time = if_else(
          is.na(time - hours(8)),
          ymd_hms(time)-hours(8),
          time-hours(8)
        )
      ) %>% 
      select(time, everything()) %>% 
      mutate(time = ymd(time)) %>% 
      group_by(time) %>% 
      summarise_all(mean, na.rm = TRUE)
  }
  
  ts_to_sum_day8a8 = function(x){
    x %>% 
      mutate(
        time = if_else(
          is.na(time - hours(8)),
          ymd_hms(time)-hours(8),
          time-hours(8)
        )
      ) %>% 
      select(time, everything()) %>% 
      mutate(time = ymd(time)) %>% 
      group_by(time) %>% 
      summarise_all(mean, na.rm = TRUE)
  }
  
  ts_to_mean_day024 = function(x){
    x %>% 
      mutate(
        time = floor_date(time, unit = "day")
      ) %>%  
      select(time, everything()) %>% 
      mutate(time = ymd(time)) %>% 
      group_by(time) %>% 
      summarise_all(mean, na.rm = TRUE)
  }
  
  ts_to_sum_day024 = function(x){
    x %>% 
      mutate(
        time = floor_date(time, unit = "day")
      ) %>%
      select(time, everything()) %>% 
      mutate(time = ymd(time)) %>% 
      group_by(time) %>% 
      summarise_all(mean, na.rm = TRUE)
  }
  
  if (fun == "mean"){
    if (mode == "8a8"){
      y = ts_to_mean_day8a8(x)
    }else if (mode == "024"){
      y = ts_to_mean_day024(x)
    }
  }else if (func == "sum"){
    if (mode == "8a8"){
      y = ts_to_sum_day8a8(x)
    }else if (mode == "024"){
      y = ts_to_sum_day024(x)
    }
  }else{
    message("select fun = mean or sum")
    return(NULL)
  }
}

calculate_error_metrics <- function(obs, sim, ...) {
  # browser()
  N <- length(obs)
  bias <- mean(sim - obs)
  pbias <- hydroGOF::pbias(sim, obs)
  rmse <- rmse(sim, obs)
  correlation <- cor(sim, obs)
  ubrmse <- sqrt(mean((sim - mean(sim) - (obs - mean(obs)))^2))
  nse <- NSE(sim, obs)
  kge_components <- KGE(sim, obs, method="2009", out.type = "full")
  kge <- kge_components$KGE.value
  alpha <- kge_components$KGE.elements[["Alpha"]]
  beta <- kge_components$KGE.elements[["Beta"]]
  r <- kge_components$KGE.elements[["r"]]
  hbe <- (sum(sim)-sum(obs))/sum(obs)
  lnNSE =   1 - sum((log(obs) - log(sim))^2, na.rm = TRUE) / sum((log(obs) - mean(log(obs), na.rm = TRUE))^2, na.rm = TRUE)
  nrmse = nrmse(sim, obs, ...)
  return(tibble(
    N = N,
    bias = bias,
    rmse = rmse,
    nrmse = nrmse,
    correlation = correlation,
    ubrmse = ubrmse,
    nse = nse,
    kge = kge,
    alpha = alpha,
    beta = beta,
    r = r,
    hbe = hbe,
    pbias = pbias,
    lnNSE = lnNSE
  ))
}

mutate_season = function(x){
  x %>% 
    mutate(m = month(time), 
           # season = get_season_v(m),
           season = case_when(
             m %in% 5:10 ~ "Accumulation",
             m %in% c(11, 12, 1:4) ~ "Release"
           ))
}
