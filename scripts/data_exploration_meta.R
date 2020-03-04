# Script to explore the data for dissertation project (META_DATA)
# Dani Gargya, daniela@gargya.de
# Feb 2020

# load libraries ----
library(tidyverse) # contains dplyr, ggplot, ...
library(maps) # for mapping the flamingo data using coordinates
library(ggthemes) # for data visualisation
library(treemap) # to create treemaps (taxa boxes)
library(wesanderson) # colour palette
library(viridis)
library("RColorBrewer")
library(treemapify) # for area graph
library(ggplot2)
library(ggpubr)

# needed?
library(broom)
library(ggalt)
library(ggrepel)
library(rgbif)
library(CoordinateCleaner)
library(gridExtra)

# load data ----
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
  mutate(duration = END_YEAR - START_YEAR)  # duration of each study

# changing All to Multiple taxa
bio_ter$TAXA <- gsub("All", "Multiple taxa", bio_ter$TAXA)

# time filtering
bio_05 <- bio_ter %>% 
  filter(START_YEAR > 1989)
# 27 from 2005, 49 from 2000

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

# spatial distribution of biodiversity time-series ----
colours1 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D")

(map_studies <- ggplot(bio_ter,
                      aes(x = CENT_LONG, y = CENT_LAT, colour = TAXA, size = NUMBER_OF_SAMPLES), alpha = I(0.7)) +
    borders("world", colour = "gray88", fill = "gray88", size = 0.3) +
    coord_cartesian(xlim = NULL, ylim = NULL, expand = TRUE) +
    theme_map() +
    geom_point(range = c(3,10)) +
    scale_colour_brewer(palette = "Dark2") +
    #scale_size(range=c(3, 10)) +
    scale_fill_manual(labels = c("Terrestrial plants",
                                 "Birds",
                                 "Mammals",
                                 "Terrestrial invertebrates",
                                 "Reptiles",
                                 "Amphibians",
                                 "Multiple Taxa")) +
    labs(title = ("\n\n a) Spatial distribution of time-series\n")) +
    theme(legend.position= "bottom", 
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.justification = "top",
          plot.title = element_text(size = 14, hjust = 0.5, face = "bold")))

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
                                                "Multiple Taxa"),
                       labels = c(1,2,3,4,5,6,7))


bio_ter$sort <- paste0(bio_ter$sort, bio_ter$START_YEAR)
bio_ter$sort <- as.numeric(as.character(bio_ter$sort))

(timeline_studies <- ggplot() +
    geom_linerange(data = bio_ter, aes(ymin = START_YEAR, ymax = END_YEAR, 
                                                  colour = TAXA,
                                                  x = fct_reorder(STUDY_ID, desc(sort))),
                   size = 1) +
    scale_colour_brewer(palette = "Dark2") +
    labs(x = NULL, y = NULL,
         title = ("\n\n b) Temporal distribution of time-series\n")) +
    theme_clean() +
    coord_flip() +
    guides(colour = F) +
    theme_clean() +
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
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 16), 
      axis.title = element_text(size = 20)))

ggsave(timeline_s, filename = "outputs/timeline_studies.png",
       height = 5, width = 8)



# taxonomic distribution of biodiversity time-series ----
# calculating sample size for each taxa
taxa_sum <- bio_ter %>%  group_by(TAXA) %>% tally

(taxa_studies <- ggplot(taxa_sum, aes(area = n, fill = TAXA, label = n,
                                   subgroup = TAXA)) +
    geom_treemap() +
    geom_treemap_subgroup_border(colour = "white", size = 1) +
    geom_treemap_text(colour = "white", place = "center", reflow = T) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(title = ("\n\n c) Taxonomic distribution of time-series\n")) +
    theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) +
    guides(fill= FALSE))

ggsave(bio_ter_area, filename = "outputs/taxa_studies.png",
       height = 5, width = 8)

# sample sizes in categories ----
# studies: multiple taxa 2, amphibians 2, terrestrial birds 35, terrestrial mammals 22, reptiles 3, terrestrial invertebrates 21, terrestrial plants 96
table_sample_size <- bio_ter %>%  
  group_by(TAXA) %>% 
  tally %>% 
  rename(Timeseries = n, Taxa = TAXA)

write.table(table_sample_size, "outputs/table_taxa.txt")

# panel ----
panel_b <- ggarrange(timeline_studies, taxa_studies, ncol = 2, align = c("h"))


panel_full <- grid.arrange(map_studies, panel_b, nrow = 2)

ggsave(panel_full, filename ="outputs/panel_studies.png",
       height = 10, width = 8)


#### draw a basic world map, add "y" or "n" for display of tropics and polar latitudes

drawWorld<-function(lats) {
  world_map<-map_data("world")
  
  g1<-ggplot()+coord_fixed()+xlab("")+ylab("")
  g1<-g1+geom_polygon(data=world_map, aes(x=long, y=lat, group=group), colour="gray60", fill="gray60")
  g1<-g1+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
               panel.background=element_rect(fill="white", colour="white"), axis.line=element_line(colour="white"),
               legend.position="none",axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())
  
  if(lats=="y") {
    g1<-g1+geom_hline(yintercept=23.5, colour="red")+geom_hline(yintercept =-23.5, colour="red")
    g1<-g1+geom_hline(yintercept=66.5, colour="darkblue")+geom_hline(yintercept =-66.5, colour="darkblue")
  }
  else { return(g1) }
  return(g1)
}

taxaCol<-c('#ffffff','#ffffbf','#5e4fa2','#f46d43','#3288bd','#abdda4','#a8c614','#d53e4f','#66c2a5','#e6f598','#fee08b','#9e0142','#fdae61')

points<-drawWorld("n")+geom_point(data=bio_ter, aes(x=CENT_LONG, y=CENT_LAT, colour=TAXA, size=TOTAL), alpha=I(0.7))
points<-points+scale_colour_manual(values=taxaCol)+scale_size(range=c(3, 10))
points
