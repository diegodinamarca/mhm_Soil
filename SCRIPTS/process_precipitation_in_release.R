pp = "G:/mHM/BASIN_DATA/data/7339001/meteo_cr2met_v2.5_1960_2023/pre.nc"
pp = rast(pp)
names(pp) = time(pp)

r.month = daily_to_monthly(pp, dates.seq = names(pp), fun = "sum")
names(r.month) = floor_date(ymd(names(pp)), unit = "months") %>% unique

month_range_mean <- function(r, range) {
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
  }

r.re = month_range_mean(r.month, 7:12)
pp.re = r.re
r.re = mask(r.re, vect(cuenca)) 
plot(r.re)
dir.create(here("DATA","RAST","PP"))
writeRaster(r.re, filename = here("DATA","RAST","PP","PP_release_mean.tif"), overwrite = TRUE)

pet = "G:/mHM/BASIN_DATA/data/7339001/meteo_cr2met_v2.5_1960_2023/pet.nc"
pet = rast(pet)
names(pet) = time(pet)

r.month = daily_to_monthly(pet, dates.seq = names(pet), fun = "sum")
names(r.month) = floor_date(ymd(names(pet)), unit = "months") %>% unique

r.re = month_range_mean(r.month, 7:12)
pet.re = r.re
plot(r.re)

mask(r.re, vect(cuenca)) %>% plot
