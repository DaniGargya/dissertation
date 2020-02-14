# Exploration data for dissertation project (FULL_DATA)
# Dani Gargya, daniela@gargya.de
# Feb 2020

# Load data ----
biotime_full <- read.csv("data/BioTIMEQuery02_04_2018.csv")

# Load libraries ----
library(tidyverse) # contains dplyr, ggplot, ...
library(maps) # for mapping the flamingo data using coordinates
library(ggthemes) # for data visualisation
library(treemap) # to create treemaps (taxa boxes)
library(treemapify) # for area graph

# data manipulation ----
bio <- biotime_full %>% 
  select(STUDY_ID, YEAR, ID_SPECIES, LATITUDE, LONGITUDE, GENUS_SPECIES)

