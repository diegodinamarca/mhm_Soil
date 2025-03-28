
# Calibracion de datos observados -----------------------------------------
datadiv = read_csv(here("SOIL MOISTURE DATA/DATA/Datos_diviner_diario_hasta2023.csv")) 
datadiv

datadiv = datadiv %>% 
  mutate(time = mdy(time)) %>% 
  rename(sitename = Nombre) %>% 
  group_by(sitename, depth) %>%
  complete(time =
             seq.Date(
               from = as.Date("2017-01-01"),
               to = as.Date("2022-12-31"),
               by = "day"
             ))
datadiv

# Eliminar depth = 10 de estos sitios porque tienen valores erroneos
# que no siguen la dinamica temporal del resto de los datos.
# puede que el tubo esten por sobre los 10cm
datafix1 = datadiv %>% 
  # rowwise() %>% 
  filter(str_detect(sitename, pattern = "BN") &
           depth != 10)

datadiv = datadiv %>% 
  filter(!str_detect(sitename, pattern = "BN")) %>% 
  bind_rows(datafix1)
datadiv

# Promediar datos diarios de diviner que tienen repeticiones
datadiv = datadiv %>% 
  group_by(sitename, depth, time) %>% 
  summarise(SM = mean(SM))
datadiv

calibrate_diviner_data <- function(x) {
  x %>% 
    ungroup() %>% 
    mutate(LC.abb = case_when(
      str_starts(sitename, pattern = "BN") ~ "NF",
      str_starts(sitename, pattern = "PI") ~ "P",
      str_starts(sitename, pattern = "MAT") ~ "S",
      str_starts(sitename, pattern = "VI") ~ "V",
      str_starts(sitename, pattern = "PR") ~ "Pr",
      str_starts(sitename, pattern = "EU") ~ "EU"
    )) %>% 
    mutate(
      SM.cal = case_when(
        LC.abb == "S" ~ SM * 1.2413 - 0.0183,
        LC.abb == "P" ~ 1.1687 * SM + 0.0876,
        LC.abb == "NF" ~ 0.0953 * log(SM) + 0.3821
      )) %>%
    filter(LC.abb == "S" |
             LC.abb == "P" |
             LC.abb == "NF")
}
datadiv_cal = calibrate_diviner_data(datadiv)
datadiv_cal %>% 
  drop_na() %>% 
  arrange(sitename, time, depth) %>%
  print(n=100)

sites_filtered_out = c("BN-MOLCO-ALTO","BN-PLAYA.BLANCA.SB","BN-SAN-ESTEBAN-1",
                         "BN-SAN-ESTEBAN-2",  "MAT-PAINOME","PI-BEBE","PINO-LA-GRANJA",
                         "PINO-MOLCO-ALTO","PINO-MOLCO-XX","PINO-SAN-ESTEBAN-1","PINO-SAN-ESTEBAN-2",
                         "BN-QUIRIHUE","MAT-BERARDO","MAT-SAN-AGUSTIN")
datadiv_cal %>% filter(!(sitename %in% sites_filtered_out)) %>% 
  filter(!is.na(SM.cal)) %>% 
  mutate(layer = if_else(depth > 100, 2, 1)) %>% 
  group_by(layer) %>% 
  summarise(n = n())


datadiv_cal_2layers = datadiv_cal %>%
  mutate(layer = if_else(depth <= 100, "L1","L2")) %>% 
  arrange(sitename, time, depth) %>% 
  group_by(sitename, time, layer) %>%
  summarise(SM = mean(SM.cal))
write_csv(datadiv_cal_2layers, here("PROC","TABLE","Obs_SM_database_diviner_calibrated_2layers.CSV"))