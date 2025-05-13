
# GitHub repo: github.com/dario-ssm/phenodev_predict

# 0. Load ----
library(terra)
library(tidyverse)
#library(lubridate)
library(here)
library(rTPC)
library(nls.multstart)
library(ggthemes)
library(chillR)
library(sf)
library(viridis)
source(here("Scripts/1-functions_phenodev_predict.R"))
library(mappestRisk)

# 1. Model fitting ------------------------------------------------------

#####  a) von Schmalensee 2023 ---------------------------------------------------

schmalensee_pieris_data <- read_delim(here("data/schmalensee2023.txt"), 
                                      delim = "\t") |> 
  filter(species == "rapae",
         life.stage == "pupa") |> 
  select(temp, dev.rate) |> 
  drop_na() 

fit_models_pieris_rapae <- mappestRisk::fit_devmodels(temp = schmalensee_pieris_data$temp,
                                                      dev_rate = schmalensee_pieris_data$dev.rate,
                                                      model_name = "all")
my_species_name <- expression(~italic("Pieris rapae"))

plot_devmodels(temp = schmalensee_pieris_data$temp,
               dev_rate = schmalensee_pieris_data$dev.rate,
               fitted_parameters = fit_models_pieris_rapae,
               species = "Pieris rapae",
               life_stage = "Pupa")+
  labs(title = my_species_name)+
  khroma::scale_color_batlow(discrete = TRUE)+
  khroma::scale_fill_batlow(discrete = TRUE)+
  theme_few()+
  theme(legend.position = "none")
ggsave(here("figures/schmalensee_tpcs.png"),
       height = 2600,
       width = 2600,
       units = "px")
  

boots_pieris_rapae <- mappestRisk::predict_curves(temp = schmalensee_pieris_data$temp,
                                                  dev_rate = schmalensee_pieris_data$dev.rate,
                                                  fitted_parameters = fit_models_pieris_rapae,
                                                  model_name_2boot = c("schoolfield", "wang", "mod_polynomial",
                                                                       "ratkowsky", "thomas","lrf", "boatman",
                                                                       "lactin2", "joehnk", "briere2", "beta", 
                                                                       "mod_weibull",                                                                       "oneill", "kamykowski", "lactin1", "pawar"),
                                                  propagate_uncertainty = TRUE,
                                                  n_boots_samples = 100) 
                                                  
mappestRisk::plot_uncertainties(bootstrap_uncertainties_tpcs = boots_pieris_rapae,
                                temp = schmalensee_pieris_data$temp,
                                dev_rate = schmalensee_pieris_data$dev.rate,
                                species = "Pieris rapae",
                                life_stage = "Pupa")
# exclude unrealistic TPC shapes or unreasonable uncertainties 
selected_models_pieris_rapae <- fit_models_pieris_rapae #|> 
 # filter(!model_name %in% c("mod_polynomial", "oneill", "mod_weibull",
  #                          "ratkowsky", "schoolfield", "wang"))
plot_devmodels(temp = schmalensee_pieris_data$temp,
               dev_rate = schmalensee_pieris_data$dev.rate,
               fitted_parameters = selected_models_pieris_rapae,
               species = "Pieris rapae",
               life_stage = "Pupa")+
  labs(title = my_species_name)+
  khroma::scale_color_batlow(discrete = TRUE)+
  khroma::scale_fill_batlow(discrete = TRUE)

write_rds(selected_models_pieris_rapae,
          here("data/selected_models_schmalensee.rds"))

#####  b) Gilbert, 1987 ---------------------------------------------------

gilbert_pieris_data <- read_delim(here("data/gilbert_pupa.csv"),
                                  delim = ";") |> 
  rename(temp = temperature,
         dev.rate = devrate,
         life_stage = ...3,
         notes = ...4,
         reference = interact) |> 
  filter(life_stage == "pupa") |> 
  dplyr::select(temp, dev.rate)
  
  
fit_models_pieris_rapae <- mappestRisk::fit_devmodels(temp = gilbert_pieris_data$temp,
                                                      dev_rate = gilbert_pieris_data$dev.rate,
                                                      model_name = "all")
my_species_name <- expression(~italic("Pieris rapae"))

plot_devmodels(temp = gilbert_pieris_data$temp,
               dev_rate = gilbert_pieris_data$dev.rate,
               fitted_parameters = fit_models_pieris_rapae,
               species = "Pieris rapae",
               life_stage = "Pupa")+
  labs(title = my_species_name)+
  khroma::scale_color_batlow(discrete = TRUE)+
  khroma::scale_fill_batlow(discrete = TRUE)+
  theme_few()+
  theme(legend.position = "none")
ggsave(here("figures/gilbert_tpcs.png"),
       width = 2600,
       height = 2600,
       units = "px")



boots_pieris_rapae <- mappestRisk::predict_curves(temp = gilbert_pieris_data$temp,
                                                  dev_rate = gilbert_pieris_data$dev.rate,
                                                  fitted_parameters = fit_models_pieris_rapae,
                                                  model_name_2boot = c("oneill", "mod_weibull", "ratkowsky",
                                                                       "lrf", "thomas", "briere2", "beta", "boatman",
                                                                       "wang", "joehnk", "lactin1", "kamykowski", "flextpc"),
                                                  propagate_uncertainty = TRUE,
                                                  n_boots_samples = 100) 

mappestRisk::plot_uncertainties(bootstrap_uncertainties_tpcs = boots_pieris_rapae,
                                temp = gilbert_pieris_data$temp,
                                dev_rate = gilbert_pieris_data$dev.rate,
                                species = "Pieris rapae",
                                life_stage = "Pupa")
# exclude unrealistic TPC shapes or unreasonable uncertainties 
selected_models_pieris_rapae_gilbert <- fit_models_pieris_rapae |> 
  filter(!model_name %in% c("mod_polynomial"))
plot_devmodels(temp = gilbert_pieris_data$temp,
               dev_rate = gilbert_pieris_data$dev.rate,
               fitted_parameters = selected_models_pieris_rapae,
               species = "Pieris rapae",
               life_stage = "Pupa")+
  labs(title = my_species_name)+
  khroma::scale_color_batlow(discrete = TRUE)+
  khroma::scale_fill_batlow(discrete = TRUE)
ggsave(here("figures/gilbert_tpcs.png"),
       width = 2600,
       height = 2600,
       units = "px")

write_rds(selected_models_pieris_rapae_gilbert,
           here("data/selected_models_gilbert.rds"))

