# Exploration data for dissertation project (FULL_DATA)
# Dani Gargya, daniela@gargya.de
# Feb 2020

### Questions
# at least 2 survey points in time # how to do that??
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
library(gridExtra)

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
  filter(!STUDY_ID == 298)  %>%  # has 147201 entries??
  select(STUDY_ID, STUDY_ID_PLOT, PLOT, NUMBER_OF_SAMPLES, TOTAL, START_YEAR, END_YEAR, duration, DATA_POINTS, TAXA, CENT_LAT, CENT_LONG, NUMBER_OF_SPECIES, studies_taxa, YEAR, ID_SPECIES, sum.allrawdata.ABUNDANCE, GENUS, SPECIES, LATITUDE, LONGITUDE)

# write and save sample
sample <- bio[1:10,]
write.csv(sample, "outputs/sample_csv.csv")

# data exploration ----
unique(bio$STUDY_ID_PLOT) # 8527 without min duration; 4929 with min duration; 152130 with study 298
unique(bio$PLOT) # does not really make sense to look at it
unique(bio$STUDY_ID) # 173 without min duration; 106 with min duration
str(bio)

# spatial scale ----
# plots per study
plot_study <- bio %>% 
  group_by(STUDY_ID) %>% 
  summarise(plots =length(unique(STUDY_ID_PLOT)))

# average/std dev plots per study
mean(plot_study$plots) # 54.49
sd(plot_study$plots) # +/- 155.03


# observations per plot
observation_plot <- bio %>% 
  group_by(STUDY_ID_PLOT) %>% 
  distinct(STUDY_ID_PLOT, NUMBER_OF_SAMPLES, TOTAL)

# average/std dev observations per plot
mean(observation_plot$NUMBER_OF_SAMPLES) # 15240
sd(observation_plot$NUMBER_OF_SAMPLES) # 39845

# temporal scale ----
# duration/data points per study ID
years_study <- bio %>% 
  group_by(STUDY_ID) %>% 
  distinct(STUDY_ID, duration, DATA_POINTS)

# average/std dev years per study ID
mean(years_study$duration) # 17.55
sd(years_study$duration) # +/- 15.79

# average/std dev data points per study ID
mean(datap_study$DATA_POINTS) # 11.08
sd(datap_study$DATA_POINTS) # 12.46


# taxa scale ----
# sample sizes taxa
samples_taxa <- bio %>% 
  group_by(TAXA) %>% 
  summarise(studies=length(unique(STUDY_ID)), # study per taxa
            plots = length(unique(STUDY_ID_PLOT)), # plots per taxa
            observations = length(NUMBER_OF_SAMPLES)) # observations per taxa
  
write.table(samples_taxa, "outputs/samples_taxa.txt")


# visualisation ----
# spatial distribution of biodiversity time-series ----
bio_short <- bio %>% 
  distinct(STUDY_ID_PLOT, CENT_LAT, CENT_LONG, TAXA, TOTAL, START_YEAR, END_YEAR, NUMBER_OF_SAMPLES)


(map_studies2 <- ggplot(bio_short,
                       aes(x = CENT_LONG, y = CENT_LAT, colour = TAXA, size = NUMBER_OF_SAMPLES), alpha = I(0.7)) +
   borders("world", colour = "gray88", fill = "gray88", size = 0.3) +
   coord_cartesian(xlim = NULL, ylim = NULL, expand = TRUE) +
   theme_map() +
   geom_point(range = c(7,15)) +
   scale_size_continuous(range = c(3,10)) +
   guides(size = FALSE) +
   scale_colour_brewer(palette = "Dark2") +
   scale_fill_manual(labels = c("Terrestrial plants",
                                "Birds",
                                "Mammals",
                                "Terrestrial invertebrates")) +
   labs(title = ("\n\n a) Spatial distribution of time-series\n")) +
   theme(legend.position= "bottom", 
         legend.title = element_blank(),
         legend.text = element_text(size = 14),
         legend.justification = "top",
         plot.title = element_text(size = 14, hjust = 0.5, face = "bold")))

ggsave(map_studies2, filename = "outputs/map_studies2.png",
       height = 5, width = 8)

# temporal distribution of biodiversity time-series ----
# making id variable as factor
bio_short$STUDY_ID_PLOT <- as.factor(as.character(bio_short$STUDY_ID_PLOT))

# create a sorting variable
bio_short$sort <- bio_short$TAXA
bio_short$sort <- factor(bio_short$sort, levels = c("Terrestrial plants",
                                                "Birds",
                                                "Mammals",
                                                "Terrestrial invertebrates"),
                       labels = c(1,2,3,4))


bio_short$sort <- paste0(bio_short$sort, bio_short$START_YEAR)
bio_short$sort <- as.numeric(as.character(bio_short$sort))

(timeline_studies2 <- ggplot() +
    geom_linerange(data = bio_short, aes(ymin = START_YEAR, ymax = END_YEAR, 
                                       colour = TAXA,
                                       x = fct_reorder(STUDY_ID_PLOT, desc(sort))),
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

ggsave(timeline_studies2, filename = "outputs/timeline_studies2.png",
       height = 5, width = 8)

# taxonomic distribution of biodiversity time-series ----
# calculating sample size for each taxa
taxa_sum <- bio_short %>%  group_by(TAXA) %>% tally

(taxa_studies2 <- ggplot(taxa_sum, aes(area = n, fill = TAXA, label = n,
                                      subgroup = TAXA)) +
    geom_treemap() +
    geom_treemap_subgroup_border(colour = "white", size = 1) +
    geom_treemap_text(colour = "white", place = "center", reflow = T) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(title = ("\n\n c) Taxonomic distribution of time-series\n")) +
    theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold")) +
    guides(fill= FALSE))

ggsave(taxa_studies2, filename = "outputs/taxa_studies2.png",
       height = 5, width = 8)

# panel ----
panel_b2 <- ggarrange(timeline_studies2, taxa_studies2, ncol = 2, align = c("h"))


panel_full2 <- grid.arrange(map_studies2, panel_b2, nrow = 2)

ggsave(panel_full2, filename ="outputs/panel_studies2.png",
       height = 10, width = 8)
