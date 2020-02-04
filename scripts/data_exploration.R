# Script to play around with the data for dissertation project
# Dani Gargya, daniela@gargya.de, Nov/2019
# Jan 2020

# load libraries ----
library(tidyverse) # contains dplyr, ggplot, ...
library(maps) # for mapping the flamingo data using coordinates
library(ggthemes) # for data visualisation
library(treemap) # to create treemaps (taxa boxes)

# load data ----
#biotime_all <- read.csv("data/BioTIMEQuery02_04_2018.csv")

biotime_all <- read.csv("data/BioTIMEMetadata_02_04_2018.csv")

# data manipulation ----
bio_ter <- biotime_all %>% 
  filter(REALM == "Terrestrial") %>% # only terrestial species
  mutate(duration = END_YEAR - START_YEAR) # duration of each study

#bio_05 <- bio_ter %>% 
  #filter(START_YEAR > 2004)

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


# visualising exploration of data ----

# spatial distribution of biodiversity time-series
(map_studies <- ggplot(bio_ter,
                      aes(x = CENT_LONG, y = CENT_LAT, colour = TAXA)) +
    borders("world", colour = "gray60", fill = "gray88", size = 0.3) +
    coord_cartesian(xlim = NULL, ylim = NULL, expand = TRUE) +
    theme_map() +
    geom_point(size = 2) +
    theme(legend.position= c(0.05, 0.2), legend.background = element_blank(),
          plot.title = element_text(size=15, hjust=0.5)) +
    labs(title="Spatial distribution studies"))

# temporal distribution of biodiversity time-series
duration_s <- bio_ter %>% 
  group_by(TAXA, duration) %>% 
  arrange(
    TAXA,
    duration,
    desc(duration)
  )

(plot_temp <- ggplot (duration_s, aes (x=x, y=y, group = TAXA, colour = TAXA)) +
    geom_segment(aes (x=START_YEAR, xend=END_YEAR, y= STUDY_ID, yend= STUDY_ID)))


#ggplot(data, aes(x=x, y=y)) +
  #geom_segment( aes(x=x, xend=x, y=0, yend=y), color="grey") +

# Taxonomic distribution of biodiversity time-series
# treemap
treemap(bio_ter,
        index="TAXA",
        vSize="STUDY_ID",
        type="index"
)
# add numbers, change title
