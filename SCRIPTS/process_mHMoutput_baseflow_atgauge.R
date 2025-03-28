
# area de cuencas
arra_area_m2 = read_sf("G:/mHM/DATA/TABLE/CAUDAL/cau_arrayan/polygon/polygon.shp") %>%
  st_transform(32718) %>%
  st_area() %>%
  as.numeric()

cau_area_m2 = read_sf(
  here("DATA","SHP","CAUQUENES.SHP")
) %>% 
  st_transform(32718) %>% 
  st_area() %>% 
  as.numeric()

df1 = read_csv(here("PROC","TABLE", str_c("mhm_output_summary_estfluv_CLSM.csv"))) %>% mutate(model = "CLSM")
df2 = read_csv(here("PROC","TABLE", str_c("mhm_output_summary_estfluv_SG.csv"))) %>% mutate(model = "SG")
df3 = read_csv(here("PROC","TABLE", str_c("mhm_output_summary_estfluv_FAO.csv"))) %>% mutate(model = "FAO")

df = bind_rows(df1, df2,df3) %>% mutate(
  across(
    where(is.numeric),
    ~ if_else(
      sitename == "cau_desem",
      .x*cau_area_m2/(1000*24*60*60),
      .x*arra_area_m2/(1000*24*60*60)
    )
  )
) 

df.month = df %>% 
  mutate(time = floor_date(time, unit = "months")) %>% 
  group_by(time, sitename, var, model) %>% 
  summarise_all(mean)

start_date = "2020-01-01"
end_date = "2023-12-01"
ggplot(df.month %>% filter(sitename == "cau_desem"))+
  geom_line(aes(x = time, y = mean, color = model), size = 1)+
  # geom_line(aes(y = qobs, color = "Qobs"), size = 1)+
  geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd, fill = model, x = time), alpha = 0.4, show.legend = F)+
  facet_wrap(.~var, scales = "free")+
  scale_x_date(limits = c(ymd(start_date),ymd(end_date)))
  

q_ltm = df.month %>% 
  mutate(m = month(time)) %>% 
  group_by(m, model, sitename, var) %>% 
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>% 
  mutate(m = factor(m, levels = c(5:12, 1:4)))

q_ltm %>% 
  filter(sitename == "cau_desem") %>% 
  ggplot(aes(x = as.numeric(m)))+
  # facet_wrap(.~var, scales = "free_y")+
  geom_line(aes(y = mean, color = model), size = 1)+
  # geom_line(aes(y = qobs, color = "Qobs"), size = 1)+
  geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd, fill = model), alpha = 0.4, show.legend = F)+
  facet_wrap(.~var, scales = "free")+
  scale_x_continuous(breaks = 1:12, labels = month.abb[c(5:12, 1:4)])+
  # scale_color_manual(name = "Model", values = model.col)+
  # scale_fill_manual(name = "Model", values = model.col)+
  labs(x = "", y = bquote("Streamflow ["*m^3*s^-1*"]"), title = "LTM monthly baseflow")
ggsave(here("RESULTADOS","FIGS","SUMMARY_BASEFLOW","LTM_STREAMFLOW.PNG"), width = 8, height = 6)


# Caudales observados en cauquenes y arrayan
qarra = read_csv(here("proc","table","Qt_Qb_Cauquenes en arrayan.csv")) %>% 
  rename(time = date, qobs = qt_m3s) %>%
  mutate(sitename = "cau_arr")
names(qarra)[3:7] = str_c("qb_",c(0.9,0.92,0.95,0.98,0.99))

qcau = read_csv(here("proc","table","Qt_Qb_Cauquenes en desembocadura.csv")) %>% 
  rename(time = date, qobs = qt_m3s) %>% 
  mutate(sitename = "cau_desem")
names(qcau)[3:7] = str_c("qb_",c(0.9,0.92,0.95,0.98,0.99))

qarra = qarra %>% 
  select(time, sitename, Qb = qb_0.9) %>% 
  mutate(model = "Obs")
qcau = qcau %>% 
  select(time, sitename, Qb = qb_0.9) %>% 
  mutate(model = "Obs")

qb_month = df %>% 
  select(time, var, model, sitename, mean) %>% 
  pivot_wider(names_from = "var", values_from = "mean") %>%
  unnest %>% 
  mutate(Qb = Qb + Qs) %>% 
  select(time, sitename, model, Qb) %>% 
  bind_rows(qarra, qcau) %>% 
  mutate(time = floor_date(time, unit = "months")) %>% 
  group_by(time, model, sitename) %>% 
  summarise(Qb = mean(Qb, na.rm = TRUE))
qb_month  
qb_month = mutate_season(qb_month)
qb_month
qb_month = qb_month %>%
  pivot_wider(names_from = model, values_from = Qb) %>% 
  pivot_longer(cols = c(CLSM,SG,FAO), names_to = "model", values_to = "Qb.sim")   %>% 
  rename(Qb.obs = Obs)

metrics = qb_month %>%
  drop_na() %>% 
  group_by(sitename, season, model) %>% 
  do(calculate_error_metrics(.$Qb.obs, .$Qb.sim))
metrics.ac = metrics %>% filter(season == "Accumulation")
metrics.re = metrics %>% filter(season == "Release")

qb_month %>% 
  # filter(alpha == "qb_0.9") %>%
  filter(season == "Accumulation") %>% 
  ggplot(aes(y = Qb.obs, x = Qb.sim))+
  geom_point(alpha = .4)+
  geom_smooth(method = "lm")+
  # ggpubr::stat_cor(aes(y = Qb.obs, x = Qb.sim, label = after_stat(rr.label)))+
  facet_wrap(
    .~sitename+model, 
    # nrow = 2, 
    # scale = "free"
  )+
  labs(x = bquote("simulated base + slow flow["*m^3*s^-1*"]"),
       y = bquote("observed base flow["*m^3*s^-1*"]"),
       caption = "Accumulation season; alpha = 0.9")+
  tune::coord_obs_pred()+
  geom_abline(color = "red")+
  geom_text(data = metrics.ac, aes(
    x = Inf, 
    y = -Inf,
    label = paste0("KGE: ",round(kge,2), "\nRMSE: ", round(rmse,2), "\nBias: ", round(bias,2))
  ), 
  hjust = 1.1, vjust = -0.5, size = 3.5, color = "blue", inherit.aes = FALSE)

qb_month %>% 
  # filter(alpha == "qb_0.9") %>%
  filter(season == "Release") %>% 
  ggplot(aes(y = Qb.obs, x = Qb.sim))+
  geom_point(alpha = .4)+
  geom_smooth(method = "lm")+
  # ggpubr::stat_cor(aes(y = Qb.obs, x = Qb.sim, label = after_stat(rr.label)))+
  facet_wrap(
    .~sitename+model, 
    # nrow = 2, 
    # scale = "free"
  )+
  labs(x = bquote("simulated base + slow flow["*m^3*s^-1*"]"),
       y = bquote("observed base flow["*m^3*s^-1*"]"),
       caption = "Release season; alpha = 0.9")+
  tune::coord_obs_pred()+
  geom_abline(color = "red")+
  geom_text(data = metrics.re, aes(
    x = Inf, 
    y = -Inf,
    label = paste0("KGE: ",round(kge,2), "\nRMSE: ", round(rmse,2), "\nBias: ", round(bias,2))
  ), 
  hjust = 1.1, vjust = -0.5, size = 3.5, color = "blue", inherit.aes = FALSE)

# Daily metrics -----------------------------------------------------------
qb_day = df %>% 
  select(time, var, model, sitename, mean) %>% 
  pivot_wider(names_from = "var", values_from = "mean") %>%
  unnest %>% 
  mutate(Qb = Qb + Qs) %>% 
  select(time, sitename, model, Qb) %>% 
  bind_rows(qarra, qcau) 
qb_day = mutate_season(qb_day)
qb_day
qb_day = qb_day %>%
  pivot_wider(names_from = model, values_from = Qb) %>% 
  pivot_longer(cols = c(CLSM,SG,FAO), names_to = "model", values_to = "Qb.sim")   %>% 
  rename(Qb.obs = Obs)

metrics = qb_day %>%
  drop_na() %>% 
  group_by(sitename, season, model) %>% 
  do(calculate_error_metrics(.$Qb.obs, .$Qb.sim))
metrics.ac = metrics %>% filter(season == "Accumulation")
metrics.re = metrics %>% filter(season == "Release")

qb_day %>% 
  # filter(alpha == "qb_0.9") %>%
  filter(season == "Accumulation") %>% 
  ggplot(aes(y = Qb.obs, x = Qb.sim))+
  geom_point(alpha = .4)+
  geom_smooth(method = "lm")+
  # ggpubr::stat_cor(aes(y = Qb.obs, x = Qb.sim, label = after_stat(rr.label)))+
  facet_wrap(
    .~sitename+model, 
    # nrow = 2, 
    # scale = "free"
  )+
  labs(x = bquote("simulated base + slow flow["*m^3*s^-1*"]"),
       y = bquote("observed base flow["*m^3*s^-1*"]"),
       caption = "Accumulation season; alpha = 0.9")+
  tune::coord_obs_pred()+
  geom_abline(color = "red")+
  geom_text(data = metrics.ac, aes(
    x = Inf, 
    y = -Inf,
    label = paste0("KGE: ",round(kge,2), "\nRMSE: ", round(rmse,2), "\nBias: ", round(bias,2))
  ), 
  hjust = 1.1, vjust = -0.5, size = 3.5, color = "blue", inherit.aes = FALSE)

qb_day %>% 
  # filter(alpha == "qb_0.9") %>%
  filter(season == "Release") %>% 
  ggplot(aes(y = Qb.obs, x = Qb.sim))+
  geom_point(alpha = .4)+
  geom_smooth(method = "lm")+
  # ggpubr::stat_cor(aes(y = Qb.obs, x = Qb.sim, label = after_stat(rr.label)))+
  facet_wrap(
    .~sitename+model, 
    # nrow = 2, 
    # scale = "free"
  )+
  labs(x = bquote("simulated base + slow flow["*m^3*s^-1*"]"),
       y = bquote("observed base flow["*m^3*s^-1*"]"),
       caption = "Release season; alpha = 0.9")+
  tune::coord_obs_pred()+
  geom_abline(color = "red")+
  geom_text(data = metrics.re, aes(
    x = Inf, 
    y = -Inf,
    label = paste0("KGE: ",round(kge,2), "\nRMSE: ", round(rmse,2), "\nBias: ", round(bias,2))
  ), 
  hjust = 1.1, vjust = -0.5, size = 3.5, color = "blue", inherit.aes = FALSE)

