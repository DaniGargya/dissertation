# Script to play around with the data for dissertation project
# Dani Gargya, daniela@gargya.de, Nov/2019
# Jan 2020

# load libraries ----
library(tidyverse) # contains dplyr, ggplot, ...
library(maps) # for mapping the flamingo data using coordinates
library(ggthemes) # for data visualisation
library(treemap) # to create treemaps (taxa boxes)
library(wesanderson) # colour palette
library(treemapify) # for area graph
library(ggplot2)

# needed?
library(broom)
library(ggalt)
library(ggrepel)
library(rgbif)
library(CoordinateCleaner)
library(gridExtra)

# load data ----
#biotime_all <- read.csv("data/BioTIMEQuery02_04_2018.csv")

biotime_all <- read.csv("data/BioTIMEMetadata_02_04_2018.csv")

# setting a clean theme ----
theme_clean <- function(){
  theme_bw() +
    theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 14, face = "plain"),             
          axis.title.y = element_text(size = 14, face = "plain"),             
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size = 15, vjust = 1, hjust = 0.5),
          legend.text = element_text(size = 12, face = "italic"),          
          legend.title = element_blank(),                              
          legend.position = c(0.5, 0.8))
}

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

# spatial distribution of biodiversity time-series ----
(map_studies <- ggplot(bio_ter,
                      aes(x = CENT_LONG, y = CENT_LAT, colour = TAXA)) +
    borders("world", colour = "gray60", fill = "gray88", size = 0.3) +
    coord_cartesian(xlim = NULL, ylim = NULL, expand = TRUE) +
    theme_map() +
    geom_point(size = 2) +
    theme(legend.position= c(0.05, 0.2), legend.background = element_blank(),
          plot.title = element_text(size=15, hjust=0.5)) +
    labs(title="Spatial distribution studies"))

ggsave(map_studies, filename = "outputs/map_studies.png",
       height = 5, width = 8)

# temporal distribution of biodiversity time-series ----
# making id variable as factor
bio_ter$STUDY_ID <- as.factor(as.character(bio_ter$STUDY_ID))

# create a sorting variable
bio_ter$sort <- bio_ter$TAXA
bio_ter$sort <- factor(bio_ter$sort, levels = c("Terrestrial plants",
                                                "Birds",
                                                "Mammals",
                                                "Terrestrial invertebrates",
                                                "Reptiles",
                                                "Amphibians",
                                                "All"),
                       labels = c(1,2,3,4,5,6,7))


bio_ter$sort <- paste0(bio_ter$sort, bio_ter$START_YEAR)
bio_ter$sort <- as.numeric(as.character(bio_ter$sort))

(timeline_s <- ggplot() +
    geom_linerange(data = bio_ter, aes(ymin = START_YEAR, ymax = END_YEAR, 
                                                  colour = TAXA,
                                                  x = fct_reorder(STUDY_ID, desc(sort))),
                   size = 1) +
    scale_colour_manual(values = wes_palette("BottleRocket1")) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    coord_flip() +
    guides(colour = F) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(),
          axis.ticks = element_blank(),
          legend.position = "bottom", 
          panel.border = element_blank(),
          legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 20, vjust = 1, hjust = 0),
          axis.text = element_text(size = 16), 
          axis.title = element_text(size = 20)))

ggsave(timeline_s, filename = "outputs/timeline_studies.png",
       height = 5, width = 8)



# Taxonomic distribution of biodiversity time-series ----
# calculating sample size for each taxa
taxa_sum <- bio_ter %>%  group_by(TAXA) %>% tally

(bio_ter_area <- ggplot(taxa_sum, aes(area = n, fill = TAXA, label = n,
                                   subgroup = TAXA)) +
    geom_treemap() +
    geom_treemap_subgroup_border(colour = "white", size = 1) +
    geom_treemap_text(colour = "white", place = "center", reflow = T) +
    scale_colour_manual(values = wes_palette("BottleRocket1")) +
    scale_fill_manual(values = wes_palette("BottleRocket1")))

ggsave(bio_ter_area, filename = "outputs/taxa_studies.png",
       height = 5, width = 8)
