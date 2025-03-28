library(here)
library(tidyverse)
library(hydroTSM)
library(hydroGOF)
library(ggforce)

source(here("Scripts","helper_functions_2.R"))

# color scales
# model.col = c("CLSM"="palevioletred3","SG"="cyan2","FAO"="gold3","Qobs"="black")


m = "FAO"
folder = here("BASIN_DATA","OUT","7339001","calibration_runs",str_c(m,"_calibrated"))
df = extract_sim_streamflow_iterations(folder)
write_csv(df, here("PROC","TABLE",str_c("simulated_streamflow_iterations_",m,".csv")))

df1 = read_csv(here("PROC","TABLE","simulated_streamflow_iterations_CLSM.csv")) %>% mutate(model = "CLSM")
df2 = read_csv(here("PROC","TABLE","simulated_streamflow_iterations_SG.csv")) %>% mutate(model = "SG")
df3 = read_csv(here("PROC","TABLE","simulated_streamflow_iterations_FAO.csv")) %>% mutate(model = "FAO")

df = bind_rows(df1,df2,df3) %>%
  mutate(
    m = month(date),
    season = case_when(
      m %in% 5:10 ~ "Accumulation",
      m %in% c(11, 12, 1:4) ~ "Release"
    ))


it_metrics = df %>%
  pivot_longer(cols = starts_with("qsim"), names_to = "iteration", values_to = "qsim") %>% 
  drop_na() %>% 
  group_by(iteration, model,season) %>% 
  do(calculate_error_metrics(.$qobs, .$qsim))
it_metrics

pmetrics = list()
pmetrics[[1]] = it_metrics %>% ggplot()+
  geom_boxplot(aes(x = model, y = r))+
  facet_wrap(.~season, nrow = 1)+
  labs(x = "", y = "Correlation (r)", title = "Daily streamflow")

pmetrics[[2]] = it_metrics %>% ggplot()+
  geom_boxplot(aes(x = model, y = bias))+
  facet_wrap(.~season, nrow = 1, scales = "free")+
  labs(x = "", y = bquote("Bias ["*m^3*s^-1*"]"), title = "Daily streamflow")

pmetrics[[3]] = it_metrics %>% ggplot()+
  geom_boxplot(aes(x = model, y = rmse))+
  facet_wrap(.~season, nrow = 1, scales = "free")+
  labs(x = "", y = bquote("RMSE ["*m^3*s^-1*"]"), title = "Daily streamflow")

pmetrics[[4]] = it_metrics %>% ggplot()+
  geom_boxplot(aes(x = model, y = kge))+
  facet_wrap(.~season, nrow = 1, scales = "free")+
  labs(x = "", y = "KGE", title = "Daily streamflow")


metrics.plot1 = ggpubr::ggarrange(
  pmetrics[[1]],pmetrics[[2]],pmetrics[[3]],pmetrics[[4]], 
  ncol = 2, nrow = 2, labels = c("a","b","c","d"))
metrics.plot1
ggsave(here("RESULTADOS","FIGS","METRICS_STREAMFLOW","daily_metrics.png"),
       width = 10, height = 8, plot = metrics.plot1)

pmetrics = list()
pmetrics[[1]] = it_metrics_month = df %>% 
  mutate(time_month = floor_date(date, unit = "month")) %>% 
  group_by(time_month, model, season) %>% 
  summarise(across(where(is.numeric), mean, na.rm =TRUE)) %>% 
  pivot_longer(cols = starts_with("qsim"), names_to = "iteration", values_to = "qsim") %>% 
  drop_na() %>% 
  group_by(iteration, model,season) %>% 
  do(calculate_error_metrics(.$qobs, .$qsim))
  
pmetrics[[1]] = it_metrics_month %>% ggplot()+
  geom_boxplot(aes(x = model, y = r))+
  facet_wrap(.~season, nrow = 1)+
  labs(x = "", y = "Correlation (r)", title = "Monthly streamflow")

pmetrics[[2]] = it_metrics_month %>% ggplot()+
  geom_boxplot(aes(x = model, y = bias))+
  facet_wrap(.~season, nrow = 1, scales = "free")+
  labs(x = "", y = bquote("Bias ["*m^3*s^-1*"]"), title = "Monthly streamflow")

pmetrics[[3]] = it_metrics_month %>% ggplot()+
  geom_boxplot(aes(x = model, y = rmse))+
  facet_wrap(.~season, nrow = 1, scales = "free")+
  labs(x = "", y = bquote("RMSE ["*m^3*s^-1*"]"), title = "Monthly streamflow")

pmetrics[[4]] = it_metrics_month %>% ggplot()+
  geom_boxplot(aes(x = model, y = kge))+
  facet_wrap(.~season, nrow = 1, scales = "free")+
  labs(x = "", y = "KGE",  title = "Monthly streamflow")
metrics.plot2 = ggpubr::ggarrange(
  pmetrics[[1]],pmetrics[[2]],pmetrics[[3]],pmetrics[[4]], 
  ncol = 2, nrow = 2, labels = c("a","b","c","d"))
metrics.plot2
ggsave(here("RESULTADOS","FIGS","METRICS_STREAMFLOW","monthly_metrics.png"), 
       width = 10, height = 8, plot = metrics.plot2)

summary_df = df %>% 
  drop_na() %>%
  filter(date > "1986-12-31",
         date < "2024-01-01") %>% 
  pivot_longer(cols = starts_with("qsim"), names_to = "iteration", values_to = "qsim") %>% 
  group_by(date, model) %>% 
  summarise(
    qobs = mean(qobs),
    qsim_min = min(qsim),
    qsim_p10 = quantile(qsim, 0.1),
    qsim_mean = mean(qsim),
    qsim_sd = sd(qsim),
    qsim_median = median(qsim),
    qsim_p90 = quantile(qsim, 0.9),
    qsim_max = max(qsim)
  )
summary_df

summary_month = summary_df %>% 
  mutate(time_month = floor_date(date, unit = "month")) %>% 
  group_by(time_month, model) %>% 
  summarise(across(where(is.numeric), mean, na.rm =TRUE))

start_date = "2020-01-01"
end_date = "2023-12-01"
summary_month %>% 
ggplot(aes(x = time_month))+
  # facet_wrap(.~var, scales = "free_y")+
  geom_line(aes(y = qsim_mean, color = model), size = 1)+
  geom_line(aes(y = qobs, color = "Qobs"), size = 1)+
  geom_ribbon(aes(ymin = qsim_mean-qsim_sd, ymax = qsim_mean+qsim_sd, fill = model), alpha = 0.4, show.legend = F)+
  scale_x_date(limits = c(ymd(start_date),ymd(end_date)))
  # scale_color_manual(name = "Model", values = model.col)+
  # scale_fill_manual(name = "Model", values = model.col)

q_ltm = summary_month %>% 
  mutate(m = month(time_month)) %>% 
  group_by(m, model) %>% 
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>% 
  mutate(m = factor(m, levels = c(5:12, 1:4)))

q_ltm %>% 
ggplot(aes(x = as.numeric(m)))+
  # facet_wrap(.~var, scales = "free_y")+
  geom_line(aes(y = qsim_mean, color = model), size = 1)+
  geom_line(aes(y = qobs, color = "Qobs"), size = 1)+
  geom_ribbon(aes(ymin = qsim_mean-qsim_sd, ymax = qsim_mean+qsim_sd, fill = model), alpha = 0.4, show.legend = F)+
  scale_x_continuous(breaks = 1:12, labels = month.abb[c(5:12, 1:4)])+
  facet_zoom(xlim = 7:12, ylim = c(0,5), horizontal = FALSE, show.area = TRUE, zoom.size = 1)+
  # scale_color_manual(name = "Model", values = model.col)+
  # scale_fill_manual(name = "Model", values = model.col)+
  labs(x = "", y = bquote("Streamflow ["*m^3*s^-1*"]"), title = "LTM monthly streamflow")
ggsave(here("RESULTADOS","FIGS","SUMMARY_STREAMFLOW","LTM_STREAMFLOW.PNG"), width = 8, height = 6)

annual_q = summary_month %>% 
  mutate(time_year = floor_date(time_month, unit = "year")) %>% 
  group_by(time_year, model) %>% 
  summarise(across(where(is.numeric), mean, na.rm=TRUE))

start_year = "1999-01-01"
end_year = "2023-12-01"
annual_q %>% 
  ggplot(aes(x = time_year))+
  geom_line(aes(y = qsim_mean, color = model), size = 1)+
  geom_line(aes(y = qobs, color = "Qobs"), size = 1)+
  geom_ribbon(aes(ymin = qsim_mean-qsim_sd, ymax = qsim_mean+qsim_sd, fill = model), alpha = 0.4, show.legend = F)+
  scale_x_date(limits = c(ymd(start_year),ymd(end_year)))+
  # scale_color_manual(name = "Model", values = model.col)+
  # scale_fill_manual(name = "Model", values = model.col)+
  labs(x = "", y = bquote("Streamflow ["*m^3*s^-1*"]"), title = "Annual streamflow")
ggsave(here("RESULTADOS","FIGS","SUMMARY_STREAMFLOW","ANNUAL_STREAMFLOW.PNG"), width = 8, height = 6)
