# Script to play around with the data for dissertation project
# Dani Gargya, daniela@gargya.de, Nov/2019
# Jan 2020

# load libraries ----
library(tidyverse) # contains dplyr, ggplot, ...

# load data ----
#biotime_all <- read.csv("data/BioTIMEQuery02_04_2018.csv")

biotime_all <- read.csv("data/BioTIMEMetadata_02_04_2018.csv")

# data manipulation ----
bio_ter <- biotime_all %>% # only terrestial species
  filter(REALM == "Terrestrial")

bio_05 <- bio_ter %>% 
  filter(START_YEAR > 2004)

# explore data ----
head(biotime_all)
tail(biotime_all)
summary(biotime_all)
str(biotime_all)

# total observation 181

# studies per taxa
bio_05 %>%
  group_by(TAXA) %>% 
  summarise(Studies=length(unique(STUDY_ID)))
# all studies: all 2, amphibians 2, birds 35, mammals 22, reptiles 3, terrestial invertebrates 21, terrestial plants 96
# 2005 onwards: amphibians 1, birds 3, mammals 5, reptiles 1, terrestrial invertebrates 5, terrestrial plants 12

# studies per climate
bio_05 %>% 
  group_by(CLIMATE) %>% 
  summarise(Studies=length(unique(STUDY_ID)))
# polar 12, polar/temperate 1, temperate 136, temperate/tropical 3, tropical 29
# temperate 21, temperate/tropical 3, tropical 3

# studies per biomemap
bio_05 %>% 
  group_by(BIOME_MAP) %>% 
  summarise(Studies=length(unique(STUDY_ID)))
