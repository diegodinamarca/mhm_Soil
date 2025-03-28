theme_grid = theme(
  panel.background = element_rect(fill = "snow", ),
  panel.grid = element_line(color = "grey", ))

plot_raster_facet <- function(r, depth, var, limits) {
  color.scale = scale_fill_distiller(type = "seq", palette = 16, direction = 1,
                                     na.value = "transparent", 
                                     name = str_c(var," ",depth),
                                     limits = limits,
                                     oob = scales::squish)
  ggplot() +
    geom_stars(data = stars::st_as_stars(r)) +  
    color.scale+
    theme_map()+
    coord_fixed()+
    facet_wrap(.~attributes)+
    theme(legend.title = element_text(size = 12))
  
}
plot_raster <- function(r, depth, var, limits) {
  color.scale = scale_fill_distiller(type = "seq", palette = 16, direction = 1,
                                     na.value = "transparent", 
                                     name = str_c(var," ",depth),
                                     limits = limits,
                                     oob = scales::squish)
  ggplot() +
    geom_stars(data = stars::st_as_stars(r)) +  
    color.scale+
    theme_map()+
    coord_fixed()+
    theme(legend.title = element_text(size = 12))+
    theme_grid
}
get_scale_lim = function(r){
  values(r) %>%
    quantile(., c(0.1, 0.9), na.rm = TRUE) %>% 
    round(., digits = 0)
}

get_raster_values <- function(r, stat, layer, season) {
  values(r) %>% 
    as_tibble() %>% 
    pivot_longer(cols = CLSM:FAO, names_to = "model", values_to = "AWC") %>% 
    mutate(stat = stat, layer = layer, season = season) %>% 
    drop_na()
}


for (i in 1:length(models)) {
  m = models[i]
  dir = here("G:/mHM/PROC/RAST/calibrated_model_output",m,"AWC")
  outdir = here("G:/mHM/PROC/RAST/calibrated_model_output",m,"AWC_summary")
  dir.create(outdir)
  files = list.files(dir, pattern = "wilting", full.names = TRUE)
  files
  r = rast(files)

  grep("soilh=5$", names(r))
  r1 = r[[seq(1,600,6)]]
  r2 = r[[seq(2,600,6)]]
  r3 = r[[seq(3,600,6)]]
  r4 = r[[seq(4,600,6)]]
  r5 = r[[seq(5,600,6)]]
  r6 = r[[seq(6,600,6)]]
  r.up = r1+r2+r3+r4+r5
  r.deep = r6
  
  r.up.pwp = r.up
  r.deep.pwp = r.deep
  
  
  calc_mean_sd <- function(r) {
    r.mean = mean(r)
    r.dif = lapply(r, function(x){
      (x-r.mean)^2
    }) %>% rast
    r.sd = sqrt(sum(r.dif)/100)
    return(list(mean = r.mean, sd = r.sd))
  }
  
  r.deep = calc_mean_sd(r.deep)
  r.up = calc_mean_sd(r.up)
  
  
  writeRaster(r.up$mean, here(outdir, str_c("PWP_up_layer_mean_",m,".tif")), overwrite = TRUE)
  writeRaster(r.up$sd, here(outdir, str_c("PWP_up_layer_sd_",m,".tif")), overwrite = TRUE)
  writeRaster(r.deep$mean, here(outdir, str_c("PWP_deep_layer_mean_",m,".tif")), overwrite = TRUE)
  writeRaster(r.deep$sd, here(outdir, str_c("PWP_deep_layer_sd_",m,".tif")), overwrite = TRUE)
  
  files = list.files(dir, pattern = "FC", full.names = TRUE)
  files
  r = rast(files)

  grep("soilh=5$", names(r))
  r1 = r[[seq(1,600,6)]]
  r2 = r[[seq(2,600,6)]]
  r3 = r[[seq(3,600,6)]]
  r4 = r[[seq(4,600,6)]]
  r5 = r[[seq(5,600,6)]]
  r6 = r[[seq(6,600,6)]]
  r.up = r1+r2+r3+r4+r5
  r.deep = r6
  
  r.up.fc = r.up
  r.deep.fc = r.deep
  

  calc_mean_sd <- function(r) {
    r.mean = mean(r)
    r.dif = lapply(r, function(x){
      (x-r.mean)^2
    }) %>% rast
    r.sd = sqrt(sum(r.dif)/100)
    return(list(mean = r.mean, sd = r.sd))
  }
  
  r.deep = calc_mean_sd(r.deep)
  r.up = calc_mean_sd(r.up)
  writeRaster(r.up$mean, here(outdir, str_c("FC_up_layer_mean_",m,".tif")), overwrite = TRUE)
  writeRaster(r.up$sd, here(outdir, str_c("FC_up_layer_sd_",m,".tif")), overwrite = TRUE)
  writeRaster(r.deep$mean, here(outdir, str_c("FC_deep_layer_mean_",m,".tif")), overwrite = TRUE)
  writeRaster(r.deep$sd, here(outdir, str_c("FC_deep_layer_sd_",m,".tif")), overwrite = TRUE)
  
  r.up.awc = r.up.fc-r.up.pwp
  r.deep.awc = r.deep.fc-r.deep.pwp
  
  r.deep = calc_mean_sd(r.deep.awc)
  r.up = calc_mean_sd(r.up.awc)
  writeRaster(r.up$mean, here(outdir, str_c("AWC_up_layer_mean_",m,".tif")), overwrite = TRUE)
  writeRaster(r.up$sd, here(outdir, str_c("AWC_up_layer_sd_",m,".tif")), overwrite = TRUE)
  writeRaster(r.deep$mean, here(outdir, str_c("AWC_deep_layer_mean_",m,".tif")), overwrite = TRUE)
  writeRaster(r.deep$sd, here(outdir, str_c("AWC_deep_layer_sd_",m,".tif")), overwrite = TRUE)
  
  r.awc = sum(r.up.awc,r.deep.awc)
  r.total = calc_mean_sd(r.awc)
  writeRaster(r.total$mean, here(outdir, str_c("AWC_total_mean_",m,".tif")), overwrite = TRUE)
  writeRaster(r.total$sd, here(outdir, str_c("AWC_total_sd_",m,".tif")), overwrite = TRUE)
  
}



for (i in 1:length(models)) {
  m = models[i]
  dir = here("G:/mHM/PROC/RAST/calibrated_model_output",m,"AWC_summary")
  outdir = here("G:/mHM/PROC/RAST/calibrated_model_output",m,"AWC_summary")
  dir.create(outdir)
  files = list.files(dir, full.names = TRUE, pattern = "AWC")
  files
  r = rast(files)
  r
}

files.clsm = list.files(here("PROC/RAST/calibrated_model_output/CLSM/AWC_summary"), pattern= "AWC", full.names = TRUE)
files.fao = list.files(here("PROC/RAST/calibrated_model_output/FAO/AWC_summary"), pattern= "AWC", full.names = TRUE)
files.sg = list.files(here("PROC/RAST/calibrated_model_output/SG/AWC_summary"), pattern= "AWC", full.names = TRUE)

awc.deep.mean = rast(c(files.clsm[1],files.sg[1],files.fao[1]))
awc.deep.sd = rast(c(files.clsm[2],files.sg[2],files.fao[2]))
awc.total.mean = rast(c(files.clsm[3],files.sg[3],files.fao[3]))
awc.total.sd = rast(c(files.clsm[4],files.sg[4],files.fao[4]))
awc.up.mean = rast(c(files.clsm[5],files.sg[5],files.fao[5]))
awc.up.sd = rast(c(files.clsm[6],files.sg[6],files.fao[6]))
names(awc.deep.mean) = c("CLSM","SG","FAO")
names(awc.deep.sd) = c("CLSM","SG","FAO")
names(awc.total.mean) = c("CLSM","SG","FAO")
names(awc.total.sd) = c("CLSM","SG","FAO")
names(awc.up.mean) = c("CLSM","SG","FAO")
names(awc.up.sd) = c("CLSM","SG","FAO")

scale.lim = get_scale_lim(awc.up.mean)
plot_raster_facet(awc.up.mean, "0-100[cm]", "AWC.mean", c(scale.lim[[1]], scale.lim[[2]]))+
  theme_grid

scale.lim = get_scale_lim(awc.up.sd)
plot_raster_facet(awc.up.sd, "0-100[cm]", "AWC.sd", c(scale.lim[[1]], scale.lim[[2]]))+
  theme_grid

scale.lim = get_scale_lim(awc.deep.mean)
plot_raster_facet(awc.deep.mean, "0-100[cm]", "AWC.mean", c(scale.lim[[1]], scale.lim[[2]]))+
  theme_grid

scale.lim = get_scale_lim(awc.deep.sd)
plot_raster_facet(awc.deep.sd, "0-100[cm]", "AWC.sd", c(scale.lim[[1]], scale.lim[[2]]))+
  theme_grid

plot.list = list()
scale.lim = get_scale_lim(awc.up.mean[[1]])
plot.list[["CLSM"]][["up"]][["mean"]] =  plot_raster(awc.up.mean[[1]], "0-100[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("CLSM - mean")
scale.lim = get_scale_lim(awc.up.mean[[2]])
plot.list[["SG"]][["up"]][["mean"]] =  plot_raster(awc.up.mean[[2]], "0-100[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("SG - mean")
scale.lim = get_scale_lim(awc.up.mean[[3]])
plot.list[["FAO"]][["up"]][["mean"]] =  plot_raster(awc.up.mean[[3]], "0-100[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("FAO - mean")

scale.lim = get_scale_lim(awc.up.sd[[1]])
plot.list[["CLSM"]][["up"]][["sd"]] =  plot_raster(awc.up.sd[[1]], "0-100[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("CLSM - sd")
scale.lim = get_scale_lim(awc.up.sd[[2]])
plot.list[["SG"]][["up"]][["sd"]] =  plot_raster(awc.up.sd[[2]], "0-100[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("SG - sd")
scale.lim = get_scale_lim(awc.up.sd[[3]])
plot.list[["FAO"]][["up"]][["sd"]] =  plot_raster(awc.up.sd[[3]], "0-100[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("FAO - sd")

scale.lim = get_scale_lim(awc.deep.mean[[1]])
plot.list[["CLSM"]][["deep"]][["mean"]] =  plot_raster(awc.deep.mean[[1]], "100-200[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("CLSM - mean")
scale.lim = get_scale_lim(awc.deep.mean[[2]])
plot.list[["SG"]][["deep"]][["mean"]] =  plot_raster(awc.deep.mean[[2]], "100-200[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("SG - mean")
scale.lim = get_scale_lim(awc.deep.mean[[3]])
plot.list[["FAO"]][["deep"]][["mean"]] =  plot_raster(awc.deep.mean[[3]], "100-200[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("FAO - mean")

scale.lim = get_scale_lim(awc.deep.sd[[1]])
plot.list[["CLSM"]][["deep"]][["sd"]] =  plot_raster(awc.deep.sd[[1]], "100-200[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("CLSM - sd")
scale.lim = get_scale_lim(awc.deep.sd[[2]])
plot.list[["SG"]][["deep"]][["sd"]] =  plot_raster(awc.deep.sd[[2]], "100-200[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("SG - sd")
scale.lim = get_scale_lim(awc.deep.sd[[3]])
plot.list[["FAO"]][["deep"]][["sd"]] =  plot_raster(awc.deep.sd[[3]], "100-200[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("FAO - sd")

ggpubr::ggarrange(plot.list[["CLSM"]][["up"]][["mean"]],
                  plot.list[["SG"]][["up"]][["mean"]], 
                  plot.list[["FAO"]][["up"]][["mean"]],
                  ncol = 3)

ggpubr::ggarrange(plot.list[["CLSM"]][["up"]][["sd"]],
                  plot.list[["SG"]][["up"]][["sd"]], 
                  plot.list[["FAO"]][["up"]][["sd"]],
                  ncol = 3)

ggpubr::ggarrange(plot.list[["CLSM"]][["deep"]][["mean"]],
                  plot.list[["SG"]][["deep"]][["mean"]], 
                  plot.list[["FAO"]][["deep"]][["mean"]],
                  ncol = 3)

ggpubr::ggarrange(plot.list[["CLSM"]][["deep"]][["sd"]],
                  plot.list[["SG"]][["deep"]][["sd"]], 
                  plot.list[["FAO"]][["deep"]][["sd"]],
                  ncol = 3)




awc.total.mean = rast(c(files.clsm[1],files.sg[1],files.fao[1]))
awc.total.sd = rast(c(files.clsm[2],files.sg[2],files.fao[2]))
awc.up.mean = rast(c(files.clsm[3],files.sg[3],files.fao[3]))
awc.up.sd = rast(c(files.clsm[4],files.sg[4],files.fao[4]))
names(awc.deep.mean) = c("CLSM","SG","FAO")
names(awc.deep.sd) = c("CLSM","SG","FAO")
names(awc.up.mean) = c("CLSM","SG","FAO")
names(awc.up.sd) = c("CLSM","SG","FAO")

scale.lim = get_scale_lim(awc.up.mean)
plot_raster_facet(awc.up.mean, "0-100[cm]", "AWC.mean", c(scale.lim[[1]], scale.lim[[2]]))+
  theme_grid

scale.lim = get_scale_lim(awc.up.sd)
plot_raster_facet(awc.up.sd, "0-100[cm]", "AWC.sd", c(scale.lim[[1]], scale.lim[[2]]))+
  theme_grid

scale.lim = get_scale_lim(awc.deep.mean)
plot_raster_facet(awc.deep.mean, "0-100[cm]", "AWC.mean", c(scale.lim[[1]], scale.lim[[2]]))+
  theme_grid

scale.lim = get_scale_lim(awc.deep.sd)
plot_raster_facet(awc.deep.sd, "0-100[cm]", "AWC.sd", c(scale.lim[[1]], scale.lim[[2]]))+
  theme_grid

plot.list = list()
scale.lim = get_scale_lim(awc.up.mean[[1]])
plot.list[["CLSM"]][["up"]][["mean"]] =  plot_raster(awc.up.mean[[1]], "0-100[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("CLSM - mean")
scale.lim = get_scale_lim(awc.up.mean[[2]])
plot.list[["SG"]][["up"]][["mean"]] =  plot_raster(awc.up.mean[[2]], "0-100[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("SG - mean")
scale.lim = get_scale_lim(awc.up.mean[[3]])
plot.list[["FAO"]][["up"]][["mean"]] =  plot_raster(awc.up.mean[[3]], "0-100[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("FAO - mean")

scale.lim = get_scale_lim(awc.up.sd[[1]])
plot.list[["CLSM"]][["up"]][["sd"]] =  plot_raster(awc.up.sd[[1]], "0-100[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("CLSM - sd")
scale.lim = get_scale_lim(awc.up.sd[[2]])
plot.list[["SG"]][["up"]][["sd"]] =  plot_raster(awc.up.sd[[2]], "0-100[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("SG - sd")
scale.lim = get_scale_lim(awc.up.sd[[3]])
plot.list[["FAO"]][["up"]][["sd"]] =  plot_raster(awc.up.sd[[3]], "0-100[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("FAO - sd")

scale.lim = get_scale_lim(awc.deep.mean[[1]])
plot.list[["CLSM"]][["deep"]][["mean"]] =  plot_raster(awc.deep.mean[[1]], "100-200[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("CLSM - mean")
scale.lim = get_scale_lim(awc.deep.mean[[2]])
plot.list[["SG"]][["deep"]][["mean"]] =  plot_raster(awc.deep.mean[[2]], "100-200[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("SG - mean")
scale.lim = get_scale_lim(awc.deep.mean[[3]])
plot.list[["FAO"]][["deep"]][["mean"]] =  plot_raster(awc.deep.mean[[3]], "100-200[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("FAO - mean")

scale.lim = get_scale_lim(awc.deep.sd[[1]])
plot.list[["CLSM"]][["deep"]][["sd"]] =  plot_raster(awc.deep.sd[[1]], "100-200[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("CLSM - sd")
scale.lim = get_scale_lim(awc.deep.sd[[2]])
plot.list[["SG"]][["deep"]][["sd"]] =  plot_raster(awc.deep.sd[[2]], "100-200[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("SG - sd")
scale.lim = get_scale_lim(awc.deep.sd[[3]])
plot.list[["FAO"]][["deep"]][["sd"]] =  plot_raster(awc.deep.sd[[3]], "100-200[cm]","AWC", c(scale.lim[[1]], scale.lim[[2]]))+ggtitle("FAO - sd")

ggpubr::ggarrange(plot.list[["CLSM"]][["up"]][["mean"]],
                  plot.list[["SG"]][["up"]][["mean"]], 
                  plot.list[["FAO"]][["up"]][["mean"]],
                  ncol = 3)

ggpubr::ggarrange(plot.list[["CLSM"]][["up"]][["sd"]],
                  plot.list[["SG"]][["up"]][["sd"]], 
                  plot.list[["FAO"]][["up"]][["sd"]],
                  ncol = 3)

ggpubr::ggarrange(plot.list[["CLSM"]][["deep"]][["mean"]],
                  plot.list[["SG"]][["deep"]][["mean"]], 
                  plot.list[["FAO"]][["deep"]][["mean"]],
                  ncol = 3)

ggpubr::ggarrange(plot.list[["CLSM"]][["deep"]][["sd"]],
                  plot.list[["SG"]][["deep"]][["sd"]], 
                  plot.list[["FAO"]][["deep"]][["sd"]],
                  ncol = 3)
