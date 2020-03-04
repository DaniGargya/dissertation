# Exploration data for dissertation project (FULL_DATA)
# Dani Gargya, daniela@gargya.de
# Feb 2020

### Questions
# at least 2 survey points in time # how to do that??
# why difference between study_id_plot and plot not bigger?
# how to get to plot and observation?
# keep 5 years min duration? omits many entries so that it does not meet min 20 studies anymore
# duration vs data points

# Load data ----
biotime_full <- read.csv("data/BioTIMEQuery02_04_2018.csv")
biotime_meta <- read.csv("data/BioTIMEMetadata_02_04_2018.csv")

# Load libraries ----
library(tidyverse) # contains dplyr, ggplot, ...
library(maps) # for mapping the flamingo data using coordinates
library(ggthemes) # for data visualisation
library(treemap) # to create treemaps (taxa boxes)
library(RColorBrewer)
library(treemapify) # for area graph
library(ggplot2)
library(ggpubr)
library(mapdata)

# data inclusion criteria and merge datasets ----

### data inclusion criteria
# terrestrial realm
# minimum time-series duration of 5 years
# at least 2 survey points in time # how to do that?? # isn't that automatic?
# at least 20 studies per taxa. 

bio <- biotime_meta %>% 
  # data inclusion criteria
  filter(REALM == "Terrestrial") %>% # only terrestial species
  group_by(STUDY_ID) %>% 
  mutate(duration = END_YEAR - START_YEAR) %>% 
  #filter(!duration < 5) %>% # minimum duration of 5 years
  #filter(!DATA_POINTS < 5) %>% 
  group_by(TAXA) %>% 
  mutate(studies_taxa=length(unique(STUDY_ID))) %>% 
  filter(!studies_taxa < 20) %>%  # minimum 20 studies per taxa
  # merge datasets
  left_join(biotime_full, by = "STUDY_ID") %>% 
  group_by(STUDY_ID) %>% 
  unite(STUDY_ID_PLOT, STUDY_ID, PLOT, sep = "_", remove=F) %>% 
  filter(!STUDY_ID == 298)  # had 147201 entries??


# data exploration ----
unique(bio$STUDY_ID_PLOT) # 8527 without min duration; 4929 with min duration; 152130 with study 298
unique(bio$PLOT) # 6888, but maybe some named the same and thats why?
unique(bio$STUDY_ID) # 173 without min duration; 106 with min duration
str(bio)


# spatial scale ----
# plots per study
plot_study <- bio %>% 
  group_by(STUDY_ID) %>% 
  summarise(plots =length(unique(STUDY_ID_PLOT)))

# average/std dev plots per study
mean(plot_study$plots) # 52.9
sqrt(var(plot_study$plots)) # +/- 152.19


# observations per plot # does not work
observation_plot <- bio %>% 
  group_by(STUDY_ID_PLOT) %>% 
  distinct(STUDY_ID_PLOT, NUMBER_OF_SAMPLES, TOTAL)

# average/std dev observations per plot
mean(observation_plot$NUMBER_OF_SAMPLES) # 15240
sqrt(var(observation_plot$NUMBER_OF_SAMPLES)) # +/- 39845.36


# temporal scale ----
# duration per study ID
years_study <- bio %>% 
  group_by(STUDY_ID) %>% 
  mutate(duration = END_YEAR - START_YEAR) %>% 
  distinct(STUDY_ID, duration)

# average/std dev years per study ID
mean(years_study$duration) # 17.21
sqrt(var(years_study$duration)) # +/- 15.61

datap_study <- bio %>%
  group_by(STUDY_ID) %>% 
  distinct(STUDY_ID, DATA_POINTS)

# average/std dev years per study ID
mean(datap_study$DATA_POINTS) # 11.08
sqrt(var(datap_study$DATA_POINTS)) # 12.46


# taxa scale ----
# study ID per taxa
study_taxa <- bio %>% 
  group_by(TAXA) %>% 
  summarise(studies=length(unique(STUDY_ID)))

# plots per taxa
plots_taxa <- bio %>% 
  group_by(TAXA) %>% 
  summarise(plots=length(unique(STUDY_ID_PLOT)))

# observations per taxa
oberservation_taxa <- bio %>% 
  group_by(TAXA) %>% 
  distinct(NUMBER_OF_SAMPLES)

# visualisation ----

