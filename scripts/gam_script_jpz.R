# GAM script for "all lakes"
# Justin Pomeranz, Fall 2023
# jpomeranz@coloradomesa.edu

# libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggeffects)
library(broom)
library(mgcv)
library(gratia)
library(broom)

# Custom functions ####
# hydr.day ####
# Custom function from Oleksy
# custom function to calculate julian-day of water year
hydro.day = function(x, start.month = 10L) {
  start.yr = year(x) - (month(x) < start.month)
  start.date = make_date(start.yr, start.month, 1L)
  as.integer(x - start.date + 1L)
}

# Gam model function ####
# Fit a GAM with specified `gam_formula`  
# Make the following objects  
  # model summaries  
  # draw_plot  
  # appraise model  
# Return a list with the following:  
  #  `lake_name`  
  # `mod` = model fit  
  # `draw_plot` = results of the `draw()` function  
  # `appraise` = results of the `apraise()` function  

gam_mod_fun <- function(data, gam_formula, lake_name){
  mod <- gam(gam_formula,
             family=Gamma(link="log"),
             data = data %>%
               filter(lake == lake_name),
             correlation = corCAR1(
               form = ~ start_year),
             method = "REML")
  mod_summary <- summary(mod)
  draw_plot <- draw(mod)
  app_plot <- appraise(mod)
  return(list(lake_name = lake_name,
              mod = mod,
              summary = mod_summary,
              draw_plot = draw_plot,
              appraise = app_plot))}


# gam summary function ####
# This function summarizes the fits from `gam_mod_fun()` above 
# and displays output in a "tidy" format

gam_summary <- function(mydata, gam_formula, group_var) {
  group_var = enquo(group_var)
  
  mydata <- mydata %>% 
    group_by(!!group_var) %>% 
    nest() 
  
  mydata %>% 
    mutate(
      model = map(data,
                  ~(gam(
                    gam_formula,
                    family=Gamma(link="log"),
                    data = .,
                    correlation = corCAR1(
                      form = ~ start_year),
                    method = "REML")))) 
}

## Non-Yellowstone lakes
non_ysl <- read_csv("scripts/other_phenology.txt")

non_ysl <- non_ysl %>%
  # select necessary columns
  select(lake, lakecode, start_year, iceOn, iceOff, orig_duration) %>% 
  # filter out lakes: 
  filter(lakecode == "JK25"| # Lake Haukivesi
           lakecode == "JK02"| # Lake Kallavesi
           lakecode == "GW369"| # Lake Kallsjön
           lakecode == "JK03"| # Näsijärvi
           lakecode == "JK05"| # Päijänne
           lakecode == "JK40"| # Pielinen
           lakecode == "NG1" # Lake Baikal
  ) %>%
  # below is jp original
  # caluclate on and off julian dates for water year
  mutate(j_on_wy = hydro.day(iceOn),
         j_off_wy = hydro.day(iceOff)) %>%
  # filter out data only from 1927 to 2022
  filter(start_year >= 1927, start_year <=2022) %>%
  mutate(ice_days = j_off_wy - j_on_wy)

# yellowstone data
ysl_ice <- read.csv("scripts/YSLoff.csv") %>%
  mutate(lake = "yellowstone")

# modify data
ysl_ice <- ysl_ice %>%
  mutate(IceOnJulian_new = case_when(IceOnJulian > 365 ~ IceOnJulian - 365,
                                     TRUE ~ IceOnJulian)) %>%
  mutate(new_iceondate = ymd(parse_date_time(paste(Year, IceOnJulian_new), orders = "yj")),
         new_iceoffdate = ymd(parse_date_time(paste(Year, IceOffJulian), orders = "yj")),
         
         IceOn_fedDOY = hydro.day(new_iceondate),
         IceOff_fedDOY = hydro.day(new_iceoffdate))#,

#Get water-year for each ice phenology indicator separately, then join back together. 
ysl_iceOn <- ysl_ice %>%
  select(new_iceondate, IceOn_fedDOY) %>%
  mutate(water_year = dataRetrieval::calcWaterYear(new_iceondate)) %>%
  filter(!new_iceondate=="1952-12-24") %>%
  arrange(water_year) %>%
  group_by(water_year) %>%
  slice_tail(n=1)

ysl_iceOn_old <- ysl_ice %>%
  select(new_iceondate, IceOn_fedDOY) %>%
  mutate(water_year = dataRetrieval::calcWaterYear(new_iceondate)) %>%
  filter(!new_iceondate=="1952-12-24")

ysl_iceOff <- ysl_ice %>%
  select(new_iceoffdate, IceOff_fedDOY, AnnualMax:SnowDepth) %>%
  mutate(water_year = dataRetrieval::calcWaterYear(new_iceoffdate))

ysl_ice_clean <- full_join(ysl_iceOn, ysl_iceOff, by="water_year") %>%
  mutate(iceDuration = IceOff_fedDOY-IceOn_fedDOY) 

# overwrite YSL with day of water year variables
ysl_ice <- ysl_ice_clean %>%
  rename(start_year = water_year,
         j_on_wy = IceOn_fedDOY,
         j_off_wy = IceOff_fedDOY,
         ice_days = iceDuration) %>%
  mutate(lake = "yellowstone") %>%
  select(lake, start_year, j_on_wy, j_off_wy, ice_days) 

# combine data sets
full_data <- bind_rows(non_ysl, ysl_ice) %>%
  mutate(start_y_c = start_year - mean(start_year, na.rm = TRUE))

# save r data structure file fo future use if needed
saveRDS(full_data,
        file = here::here("scripts/all_lake_data.rds"))

# for-loop fitting models ####
lake_names <- full_data %>%
  pull(lake) %>%
  unique()

# ice-on fits ####
ice_on_list <- list()
for (lake in 1:length(lake_names)){
  lake_name = lake_names[lake]
  response_mod <- gam_mod_fun(full_data, j_on_wy~s(start_year), lake_name)
  ice_on_list[[lake]] <- response_mod
}
names(ice_on_list) <- lake_names

# ice-off fits
ice_off_list <- list()
for (lake in 1:length(lake_names)){
  lake_name = lake_names[lake]
  response_mod <- gam_mod_fun(full_data, j_off_wy~s(start_year), lake_name)
  ice_off_list[[lake]] <- response_mod
}
names(ice_off_list) <- lake_names

# ice duration fits
ice_days_list <- list()
for (lake in 1:length(lake_names)){
  lake_name = lake_names[lake]
  response_mod <- gam_mod_fun(full_data, ice_days~s(start_year), lake_name)
  ice_days_list[[lake]] <- response_mod
}
names(ice_days_list) <- lake_names

# gam summaries ####
# ice-on summary ####
on_stats <- gam_summary(full_data, j_on_wy~s(start_year), lake) %>%
  mutate(tidied = map(model, tidy)) %>%
  select(tidied) %>%
  unnest(tidied) %>%
  mutate(phenology = "ice_on") %>%
  arrange(lake_name, p.value)
on_stats

# ice-off summary ####
off_stats <- gam_summary(full_data, j_off_wy~s(start_year), lake) %>%
  mutate(tidied = map(model, tidy)) %>%
  select(tidied) %>%
  unnest(tidied) %>%
  arrange(lake_name, p.value) %>%
  mutate(phenology = "ice_off") 
off_stats

# ice duration summary ####
duration_stats <- gam_summary(full_data, ice_days~s(start_year), lake) %>%
  mutate(tidied = map(model, tidy)) %>%
  select(tidied) %>%
  unnest(tidied) %>%
  arrange(lake_name, p.value) %>%
  mutate(phenology = "ice_duration") 
duration_stats

# write csv of gam summaries ####
bind_rows(on_stats, off_stats, duration_stats) %>%
  write_csv("scripts/gam_stat_table.csv")

# Full figure ####
full_data %>% 
  select(lake, start_year, ice_days, j_on_wy, j_off_wy) %>%
  pivot_longer(ice_days:j_off_wy) %>%
  ggplot(aes(x = start_year,
             y = value)) +
  geom_point() +
  theme_bw() +
  facet_grid(name ~ lake, scales = "free_y") +
  geom_smooth(method = "gam") + 
  labs(y = "Day of water year") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("scripts/full_fig_gams.png",
       units = "in",
       height = 8,
       width = 8,
       dpi = 350)
