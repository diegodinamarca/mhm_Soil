
# Calibracion de datos observados -----------------------------------------
datasondas = read_csv(here("SOIL MOISTURE DATA/DATA/sondas_humedad_suelo_diario_day8a8.csv")) 
datasondas = datasondas %>% 
  mutate(time = mdy(time))
datasondas %>% 
  group_by(sitename) %>% 
  summarise(start = min(time),
            end = max(time),
            n = n())

# calibrated probe sites
cal_sites = c("PINO-SANPEDRO","PINO-LAGRANJA","MAT-LAGRANJA",
              "MAT-SASB","PINO-SASB","PINO-GRANIER","BN-METAMORFICO",
              "BN-PLAYABLANCA","MAT-INIA")

# calibracion de sitios
cal_x2 = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
cal_x = c(1.27, 1.21, 1.16, 1.05, 1.29, 1, 1, 1, 1.194)
cal_int = c(0.004,-0.008,-0.027,-0.008,-0.032, 0, 0, 0, -0.0162)

# write_csv(tibble(cal_sites, "ax2" = cal_x2, "bx" = cal_x, "c" = cal_int),
#           here("PROC","PROBE_CALIBRATION_EQUATIONS.CSV"))

datasondas_cal = datasondas %>% 
  filter(sitename %in% cal_sites) %>%
  mutate(SM = case_when(
    sitename == cal_sites[1] ~ (SM)^2*cal_x2[1]+SM*cal_x[1]+cal_int[1],
    sitename == cal_sites[2] ~ (SM)^2*cal_x2[2]+SM*cal_x[2]+cal_int[2],
    sitename == cal_sites[3] ~ (SM)^2*cal_x2[3]+SM*cal_x[3]+cal_int[3],
    sitename == cal_sites[4] ~ (SM)^2*cal_x2[4]+SM*cal_x[4]+cal_int[4],
    sitename == cal_sites[5] ~ (SM)^2*cal_x2[5]+SM*cal_x[5]+cal_int[5],
    sitename == cal_sites[6] ~ (SM)^2*cal_x2[6]+SM*cal_x[6]+cal_int[6],
    sitename == cal_sites[7] ~ (SM)^2*cal_x2[7]+SM*cal_x[7]+cal_int[7],
    sitename == cal_sites[8] ~ (SM)^2*cal_x2[8]+SM*cal_x[8]+cal_int[8],
    sitename == cal_sites[9] ~ (SM)^2*cal_x2[9]+SM*cal_x[9]+cal_int[9]
  ))

datasondas%>% 
  filter(sitename %in% cal_sites) %>%
  ggplot(aes(x = SM)) +
  geom_density(aes(color = "pre-calib"))+
  geom_density(data = datasondas_cal, aes(x = SM, color = "calib"))+
  facet_wrap(.~sitename, scales = "free")
# datasondas_cal %>% filter(sitename == "BN-METAMORFICO")

datasondas_cal %>% 
  mutate(layer = if_else(depth > 100, 2, 1)) %>% 
  group_by(layer) %>% 
  summarise(n = n())

####################################################################.
########### Homogenizar datos observados a horizontes mHM ###########
####################################################################.

df1_1 = datasondas_cal %>%
  filter(sitename == "PINO-SASB") %>% 
  pivot_wider(names_from = depth, values_from = SM, names_prefix = "D") %>% 
  filter(!is.na(D28)) %>% 
  mutate(L1 = NA, 
         L2 = NA, 
         L3 = D28, 
         L4 = (D58*7+D28*3)/10,
         L5 = D88,
         L6 = (D128+D178)/2) %>% 
  select(time, sitename, L1:L6)

df1_2 = datasondas_cal %>%
  filter(sitename == "PINO-SASB") %>% 
  pivot_wider(names_from = depth, values_from = SM, names_prefix = "D") %>% 
  filter(is.na(D28)) %>% 
  mutate(L1 = NA, 
         L2 = NA, 
         L3 = D30, 
         L4 = (D30+D60)/2,
         L5 = D88,
         L6 = (D100+2*D140)/3) %>% 
  select(time, sitename, L1:L6)

df2 = datasondas_cal %>%
  filter(sitename == "PINO-LAGRANJA") %>% 
  pivot_wider(names_from = depth, values_from = SM, names_prefix = "D") %>% 
  mutate(L1 = NA, 
         L2 = NA, 
         L3 = D30, 
         L4 = (D30+D60)/2,
         L5 = D90,
         L6 = (D120+D160)/2) %>% 
  select(time, sitename, L1:L6)

df3 = datasondas_cal %>%
  filter(sitename == "MAT-LAGRANJA") %>% 
  pivot_wider(names_from = depth, values_from = SM, names_prefix = "D") %>% 
  mutate(L1 = NA, 
         L2 = NA, 
         L3 = D30, 
         L4 = (D30+D60)/2,
         L5 = (2*D80+D100)/3,
         L6 = (D100+2*D170)/3) %>% 
  select(time, sitename, L1:L6)

df4_1 = datasondas_cal %>%
  filter(sitename == "MAT-SASB") %>% 
  pivot_wider(names_from = depth, values_from = SM, names_prefix = "D") %>%
  filter(is.na(D47)) %>% 
  
  mutate(L1 = NA, 
         L2 = NA, 
         L3 = NA, 
         L4 = D40,
         L5 = (D70+D90)/2,
         L6 = (D120+D170)/2) %>% 
  select(time, sitename, L1:L6)

df4_2 = datasondas_cal %>%
  filter(sitename == "MAT-SASB") %>% 
  pivot_wider(names_from = depth, values_from = SM, names_prefix = "D") %>%
  filter(!is.na(D47)) %>% 
  mutate(L1 = NA, 
         L2 = NA, 
         L3 = NA, 
         L4 = D47,
         L5 = (D77*10+D97*8)/18,
         L6 = (D97*2+D127*10+D177*10)/22) %>% 
  select(time, sitename, L1:L6)

df5 = datasondas_cal %>%
  filter(sitename == "PINO-SANPEDRO") %>% 
  pivot_wider(names_from = depth, values_from = SM, names_prefix = "D") %>% 
  mutate(L1 = NA, 
         L2 = D10, 
         L3 = D30, 
         L4 = D50,
         L5 = D90,
         L6 = D160) %>% 
  select(time, sitename, L1:L6)

df6 = datasondas_cal %>%
  filter(sitename == "PINO-GRANIER") %>% 
  pivot_wider(names_from = depth, values_from = SM, names_prefix = "D") %>% 
  mutate(L1 = NA, 
         L2 = NA, 
         L3 = NA,
         L4 = (D35*2+D60)/3,
         L5 = (D60+D100)/2,
         L6 = (D140+D180)/2) %>% 
  select(time, sitename, L1:L6)

df7 = datasondas_cal %>%
  filter(sitename == "BN-METAMORFICO") %>% 
  pivot_wider(names_from = depth, values_from = SM, names_prefix = "D") %>% 
  mutate(L1 = NA, 
         L2 = D10, 
         L3 = NA, 
         L4 = D37,
         L5 = (D67+D83)/2,
         L6 = NA) %>% 
  select(time, sitename, L1:L6)


df8 = datasondas_cal %>%
  filter(sitename == "BN-PLAYABLANCA") %>% 
  pivot_wider(names_from = depth, values_from = SM, names_prefix = "D") %>% 
  mutate(L1 = NA, 
         L2 = D10, 
         L3 = NA, 
         L4 = D40,
         L5 = D70,
         L6 = D120) %>% 
  select(time, sitename, L1:L6)

df9 = datasondas_cal %>%
  filter(sitename == "MAT-INIA") %>% 
  pivot_wider(names_from = depth, values_from = SM, names_prefix = "D") %>% 
  mutate(L1 = NA, 
         L2 = NA, 
         L3 = D20, 
         L4 = D50,
         L5 = D80,
         L6 = (D120+D150)/2) %>% 
  select(time, sitename, L1:L6)

datasondas_homo = bind_rows(df1_1, df1_2, df2, df3, 
                            df4_1, df4_2, df5, df6,
                            df7, df8, df9) %>% 
  pivot_longer(cols = L1:L6, names_to = "layer", values_to = "swc_obs") 
write_csv(datasondas_homo, here("PROC","TABLE","Obs_SM_database_probes_calibrated.CSV"))

ggplot(data = datasondas_homo, aes(x = swc_obs, color = layer))+
  geom_density(size = .8, alpha = .5)+
  # geom_density(data = datasondas_cal)
  facet_wrap(.~sitename, scales = "free_y")+
  scale_color_brewer(type = "qual", palette = 8)+
  theme_bw()+
  labs(x = "Soil Moisture (v/v)", title = "Homogenized Observed Soil Moisture",
       color = "Layer")


depths.width = c(5,10,15,30,40,100)
depths.pond = depths.width/100
names(depths.pond) = str_c("L",1:6)
datasondas_2layers = datasondas_homo %>% 
  mutate(swc_obs = depths.pond[layer]*swc_obs) %>% 
  mutate(layer = if_else(layer %in% str_c("L",1:5),1,2)) %>% 
  group_by(sitename,time, layer) %>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
  filter(swc_obs != 0)

datasondas_2layers
write_csv(datasondas_2layers, here("PROC","TABLE","Obs_SM_database_probes_calibrated_2layers.CSV"))


# ggsave(filename = here("Resultados","figs","homogenized_obs_soil_moist.png"),
#        width = 8, height = 4)