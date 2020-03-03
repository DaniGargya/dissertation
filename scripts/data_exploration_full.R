# Exploration data for dissertation project (FULL_DATA)
# Dani Gargya, daniela@gargya.de
# Feb 2020

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

# merge datasets ----
bio <- biotime_meta %>% 
  filter(REALM == "Terrestrial") %>% # only terrestial species
  left_join(biotime_full, by = "STUDY_ID") %>% 
  unite(STUDY_ID_PLOT, STUDY_ID, PLOT, sep = "_", remove=F) %>% 
  group_by(STUDY_ID) %>% 
  filter(!STUDY_ID == 298) # had 147201 entries??
  
unique(bio$STUDY_ID_PLOT) # 156 000?
unique(bio$PLOT)


# data manipulation ----
bio$TAXA <- gsub("All", "Multiple taxa", bio$TAXA)

# data exploration ----
str(bio)
unique(bio$PLOT) # 170 224??
unique (bio$STUDY_ID) # 181

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
  summarise(observations =length(unique(PLOT)))

# average/std dev observations per plot
mean(observation_plot$observations) # 52.9
sqrt(var(plot_study$plots)) # +/- 152.19


# temporal scale ----
# duration per study ID
years_study <- bio %>% 
  group_by(STUDY_ID) %>% 
  mutate(duration = END_YEAR - START_YEAR) %>% 
  distinct(STUDY_ID, duration)

# average/std dev years per study ID
mean(years_study$duration) # 17.21
sqrt(var(years_study$duration)) # +/- 15.61

# calculate duration of plots?

# taxa scale ----
# study ID per taxa
study_taxa <- bio %>% 
  group_by(TAXA) %>% 
  summarise(studies=length(unique(STUDY_ID)))

# plots per taxa
plots_taxa <- bio %>% 
  group_by(TAXA) %>% 
  summarise(plots=length(unique(STUDY_ID_PLOT)))

# observations per taxa # something is wrong
oberservations_taxa <- bio %>% 
  group_by(TAXA) %>% 
  summarise(obs=length(unique(PLOT)))
