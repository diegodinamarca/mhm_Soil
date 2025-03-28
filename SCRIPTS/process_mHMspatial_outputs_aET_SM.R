m = "CLSM"

models = c("CLSM","SG","FAO")
for (i in 2:length(models)){
  m = models[i]
  files = list.files(str_c("G:/mHM/PROC/RAST/calibrated_model_output/",m,"/aET_0.100"), full.names = TRUE)
  files
  
  # release months 7:12
  month_range_mean <- function(files, range) {
    lapply(files, function(f){
      r = rast(f)
      # resta 4 meses para comenzar el año hidrologico en mayo
      new_times = ymd(names(r)) %m-% months(4)
      new_times
      years = year(new_times) %>% unique
      years = years[2:(length(years)-1)]
      
      it.mean = lapply(years, function(y){
        index = grep(y, new_times)
        y.r = r[[index]]
        y.r = sum(y.r[[range]])
      }) %>% rast %>% mean
    })
  }
  
  r.it = month_range_mean(files, 7:12)
  r.ltm = rast(r.it)
  r.mean = mean(r.ltm)
  writeRaster(r.mean, here("PROC","RAST","CALIBRATED_MODEL_OUTPUT",m,"season_summary_ok","release_0_100_mean.tif"))
  r.sd = sqrt(sum((r.mean - r.ltm)^2)/100)
  writeRaster(r.sd, here("PROC","RAST","CALIBRATED_MODEL_OUTPUT",m,"season_summary_ok","release_0_100_sd.tif"))
  plot(r.sd)
  plot(r.mean)
  
  
  r.it = month_range_mean(files, 1:6)
  r.ltm = rast(r.it)
  r.mean = mean(r.ltm)
  writeRaster(r.mean, here("PROC","RAST","CALIBRATED_MODEL_OUTPUT",m,"season_summary_ok","accum_0_100_mean.tif"))
  r.sd = sqrt(sum((r.mean - r.ltm)^2)/100)
  writeRaster(r.sd, here("PROC","RAST","CALIBRATED_MODEL_OUTPUT",m,"season_summary_ok","accum_0_100_sd.tif"))
  plot(r.sd)
  plot(r.mean)
  
  
  # fies for the deeper soil layers
  files = list.files(str_c("G:/mHM/PROC/RAST/calibrated_model_output/",m,"/aET"), full.names = TRUE)
  index = grep("H6", files)
  files = files[index]
  files
  
  r.it = month_range_mean(files, 7:12)
  r.ltm = rast(r.it)
  r.mean = mean(r.ltm)
  writeRaster(r.mean, here("PROC","RAST","CALIBRATED_MODEL_OUTPUT",m,"season_summary_ok","release_100_200_mean.tif"))
  r.sd = sqrt(sum((r.mean - r.ltm)^2)/100)
  writeRaster(r.sd, here("PROC","RAST","CALIBRATED_MODEL_OUTPUT",m,"season_summary_ok","release_100_200_sd.tif"))
  plot(r.sd)
  plot(r.mean)
  
  
  r.it = month_range_mean(files, 1:6)
  r.ltm = rast(r.it)
  r.mean = mean(r.ltm)
  writeRaster(r.mean, here("PROC","RAST","CALIBRATED_MODEL_OUTPUT",m,"season_summary_ok","accum_100_200_mean.tif"))
  r.sd = sqrt(sum((r.mean - r.ltm)^2)/100)
  writeRaster(r.sd, here("PROC","RAST","CALIBRATED_MODEL_OUTPUT",m,"season_summary_ok","accum_100_200_sd.tif"))
  plot(r.sd)
  plot(r.mean)
}



# SWC ---------------------------------------------------------------------


for (i in 1:length(models)){
  m = models[i]
  files = list.files(str_c("G:/mHM/PROC/RAST/calibrated_model_output/",m,"/SWC_0.100"), full.names = TRUE)
  files
  
  # release months 7:12
  month_range_mean <- function(files, range) {
    lapply(files, function(f){
      r = rast(f)
      # resta 4 meses para comenzar el año hidrologico en mayo
      new_times = ymd(names(r)) %m-% months(4)
      new_times
      years = year(new_times) %>% unique
      years = years[2:(length(years)-1)]
      
      it.mean = lapply(years, function(y){
        index = grep(y, new_times)
        y.r = r[[index]]
        y.r = mean(y.r[[range]])
      }) %>% rast %>% mean
    })
  }
  
  r.it = month_range_mean(files, 7:12)
  r.ltm = rast(r.it)
  r.mean = mean(r.ltm)
  writeRaster(r.mean, here("PROC","RAST","CALIBRATED_MODEL_OUTPUT",m,"season_summary_sm","release_0_100_mean.tif"))
  r.sd = sqrt(sum((r.mean - r.ltm)^2)/100)
  writeRaster(r.sd, here("PROC","RAST","CALIBRATED_MODEL_OUTPUT",m,"season_summary_sm","release_0_100_sd.tif"))
  plot(r.sd)
  plot(r.mean)
  
  
  r.it = month_range_mean(files, 1:6)
  r.ltm = rast(r.it)
  r.mean = mean(r.ltm)
  writeRaster(r.mean, here("PROC","RAST","CALIBRATED_MODEL_OUTPUT",m,"season_summary_sm","accum_0_100_mean.tif"))
  r.sd = sqrt(sum((r.mean - r.ltm)^2)/100)
  writeRaster(r.sd, here("PROC","RAST","CALIBRATED_MODEL_OUTPUT",m,"season_summary_sm","accum_0_100_sd.tif"))
  plot(r.sd)
  plot(r.mean)
  
  
  # fies for the deeper soil layers
  files = list.files(str_c("G:/mHM/PROC/RAST/calibrated_model_output/",m,"/SM"), full.names = TRUE)
  index = grep("H6", files)
  files = files[index]
  files
  month_range_mean <- function(files, range) {
    lapply(files, function(f){
      r = rast(f)*1000
      # resta 4 meses para comenzar el año hidrologico en mayo
      new_times = ymd(names(r)) %m-% months(4)
      new_times
      years = year(new_times) %>% unique
      years = years[2:(length(years)-1)]
      
      it.mean = lapply(years, function(y){
        index = grep(y, new_times)
        y.r = r[[index]]
        y.r = mean(y.r[[range]])
      }) %>% rast %>% mean
    })
  }
  
  r.it = month_range_mean(files, 7:12)
  r.ltm = rast(r.it)
  r.mean = mean(r.ltm)
  writeRaster(r.mean, here("PROC","RAST","CALIBRATED_MODEL_OUTPUT",m,"season_summary_sm","release_100_200_mean.tif"))
  r.sd = sqrt(sum((r.mean - r.ltm)^2)/100)
  writeRaster(r.sd, here("PROC","RAST","CALIBRATED_MODEL_OUTPUT",m,"season_summary_sm","release_100_200_sd.tif"))
  plot(r.sd)
  plot(r.mean)
  
  
  r.it = month_range_mean(files, 1:6)
  r.ltm = rast(r.it)
  r.mean = mean(r.ltm)
  writeRaster(r.mean, here("PROC","RAST","CALIBRATED_MODEL_OUTPUT",m,"season_summary_sm","accum_100_200_mean.tif"))
  r.sd = sqrt(sum((r.mean - r.ltm)^2)/100)
  writeRaster(r.sd, here("PROC","RAST","CALIBRATED_MODEL_OUTPUT",m,"season_summary_sm","accum_100_200_sd.tif"))
  plot(r.sd)
  plot(r.mean)
}
