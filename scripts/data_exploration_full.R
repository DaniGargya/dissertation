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
library("RColorBrewer")
library(treemapify) # for area graph
library(ggplot2)
library(ggpubr)

# merge datasets ----
bio <- biotime_meta %>% 
  left_join(biotime_full, by = "STUDY_ID") %>% 
  filter(REALM == "Terrestrial") %>% # only terrestial species
  mutate(duration = END_YEAR - START_YEAR)

# data manipulation ----
bio$TAXA <- gsub("All", "Multiple taxa", bio$TAXA)

# data exploration ----
str(bio)
unique(bio$PLOT) # 170 224??
unique (bio$STUDY_ID) # 181