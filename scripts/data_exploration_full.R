# Exploration data for dissertation project (FULL_DATA)
# Dani Gargya, daniela@gargya.de
# Feb 2020
display.brewer.pal(n = 8, name = 'Dark2')
# Hexadecimal color specification 
brewer.pal(n = 8, name = "Dark2")

### Questions
# min duration vs data points in time-series
# exclude 502 as well? has 1451 studies
# info: study 327 has 171843 number of samples in chile

# Load data ----
biotime_full <- read.csv("data/BioTIMEQuery02_04_2018.csv")
biotime_meta <- read.csv("data/BioTIMEMetadata_02_04_2018.csv")
bio <- read.csv("data/bio.csv")

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
# minimum study points 5 years
# at least 2 survey points in time # how to do that?? # isn't that automatic?
# at least 15 studies per taxa
# no more than 5000 plots per study

bio <- biotime_meta %>% 
  # data inclusion criteria
  filter(REALM == "Terrestrial") %>% # only terrestial species
  group_by(STUDY_ID) %>% 
  mutate(duration = END_YEAR - START_YEAR) %>% 
  filter(!duration < 5) %>% # minimum duration of 5 years
  #filter(!DATA_POINTS < 5) %>% 
  group_by(TAXA) %>% 
  mutate(studies_taxa=length(unique(STUDY_ID))) %>% 
  filter(!studies_taxa < 15) %>%  # minimum 20 studies per taxa
  ungroup() %>% 
  # merge datasets
  left_join(biotime_full, by = "STUDY_ID") %>% 
  group_by(STUDY_ID) %>% 
  unite(STUDY_ID_PLOT, STUDY_ID, PLOT, sep = "_", remove=F) %>% 
  filter(!STUDY_ID == 298)  %>%  # has 147201 entries?? set upper limit
  group_by(STUDY_ID_PLOT) %>% 
  filter(HAS_PLOT == "Y") %>% 
  filter(YEAR %in% c(max(YEAR), min(YEAR))) %>% 
  mutate(number_plots = length(unique(YEAR))) %>% 
  filter(number_plots == 2) %>% 
  dplyr::select(STUDY_ID, STUDY_ID_PLOT, PLOT, SAMPLE_DESC, NUMBER_OF_SAMPLES, TOTAL, START_YEAR, END_YEAR, duration, DATA_POINTS, TAXA, CENT_LAT, CENT_LONG, LATITUDE, LONGITUDE, NUMBER_OF_SPECIES, studies_taxa, YEAR, ID_SPECIES, sum.allrawdata.ABUNDANCE, GENUS_SPECIES, LATITUDE, LONGITUDE, AREA_SQ_KM, HAS_PLOT, NUMBER_LAT_LONG) %>% 
  ungroup()

bio_short2 <- bio %>% 
  distinct(STUDY_ID_PLOT, .keep_all = TRUE)

# saving data subset
#write.csv(bio, "data/bio.csv")

# meta filtered
bio_meta <- biotime_meta %>% 
  # data inclusion criteria
  filter(REALM == "Terrestrial") %>% # only terrestial species
  group_by(STUDY_ID) %>% 
  mutate(duration = END_YEAR - START_YEAR) %>% 
  filter(!duration < 5) %>% # minimum duration of 5 years
  #filter(!DATA_POINTS < 5) %>% 
  group_by(TAXA) %>% 
  mutate(studies_taxa=length(unique(STUDY_ID))) %>% 
  filter(!studies_taxa < 15) %>% 
  filter(!STUDY_ID == 298)
  

# write and save sample
#sample <- bio[1:10,]
#write.csv(sample, "outputs/sample_csv.csv")

# data exploration ----
unique(bio$STUDY_ID_PLOT) # 7473
unique(bio$PLOT) # does not really make sense to look at it
unique(bio$STUDY_ID) # 139
str(bio)

# spatial scale ----
# plots per study
plot_study <- data1 %>% 
  group_by(STUDY_ID) %>% 
  summarise(plots =length(unique(STUDY_ID_PLOT)))

# average/std dev plots per study
mean(plot_study$plots) # 64
sd(plot_study$plots) # +/- 145


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
mean(years_study$DATA_POINTS) # 11.08
sd(years_study$DATA_POINTS) # 12.46

min(data1$START_YEAR)
max(data1$END_YEAR)
mean(data1$duration_plot)
sd(data1$duration_plot)


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

  
#distinct(STUDY_ID_PLOT, STUDY_ID, TAXA, TOTAL, START_YEAR, END_YEAR, NUMBER_OF_SAMPLES, AREA_SQ_KM, HAS_PLOT, NUMBER_LAT_LONG)

(map_studies2 <- ggplot(bio_short2,
                       aes(x = LONGITUDE, y = LATITUDE, colour = TAXA, size = NUMBER_OF_SAMPLES), alpha = I(0.7)) +
   borders("world", colour = "gray88", fill = "gray88", size = 0.3) +
   coord_cartesian(xlim = NULL, ylim = NULL, expand = TRUE) +
   theme_map() +
   geom_point(range = c(7,15)) +
   scale_size_continuous(range = c(3,10)) +
   guides(size = FALSE) +
   scale_colour_manual(values = taxa.palette) +
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
bio_short2$STUDY_ID_PLOT <- as.factor(as.character(bio_short2$STUDY_ID_PLOT))

# create a sorting variable
bio_short2$sort <- bio_short2$TAXA
bio_short2$sort <- factor(bio_short2$sort, levels = c("Terrestrial plants",
                                                "Birds",
                                                "Mammals",
                                                "Terrestrial invertebrates"),
                       labels = c(1,2,3,4))


bio_short2$sort <- paste0(bio_short2$sort, bio_short2$START_YEAR)
bio_short2$sort <- as.numeric(as.character(bio_short2$sort))

(timeline_studies2 <- ggplot() +
    geom_linerange(data = bio_short2, aes(ymin = START_YEAR, ymax = END_YEAR, 
                                       colour = TAXA,
                                       x = fct_reorder(STUDY_ID_PLOT, desc(sort))),
                   size = 1) +
    scale_colour_manual(values = taxa.palette) +
    labs(x = NULL, y = NULL,
         title = ("\n\n b) Temporal distribution of time-series\n")) +
    #theme_clean() +
    coord_flip() +
    guides(colour = F) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom", 
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
taxa_sum <- bio_short2 %>%  group_by(TAXA) %>% tally

(taxa_studies2 <- ggplot(taxa_sum, aes(area = n, fill = TAXA, label = n,
                                      subgroup = TAXA)) +
    geom_treemap() +
    geom_treemap_subgroup_border(colour = "white", size = 1) +
    geom_treemap_text(colour = "white", place = "center", reflow = T) +
    scale_colour_manual(values = taxa.palette) +
    scale_fill_manual(values = taxa.palette) +
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

# checking for need rarefaction ----
# histograms area ----
# histogram area studyID_plot
png("outputs/area_hist_base.png")
area_hist_base <- hist(bio_short$AREA_SQ_KM, breaks = 100)
dev.off()

# log axis
log_area <- log(bio_short$AREA_SQ_KM)

png("outputs/hist_log_area.png")
hist_log_area <- hist(log_area, breaks =100)
dev.off()

# histogram km2 area plots
#(area_hist3 <- ggplot(bio, aes(x = AREA_SQ_KM)) +
   #geom_histogram()+
   #theme_clean())

# histogram km2 area study ID
#(area_hist2 <- ggplot(bio_meta, aes(x = AREA_SQ_KM)) +
    #geom_histogram() +
    #theme_clean())

# histogram km2 area studyID_plot
#(area_hist <- ggplot(bio_short, aes(x = AREA_SQ_KM)) +
    #geom_histogram()+
    #theme_clean())
#ggsave(area_hist, filename = "outputs/area_hist.png", height =7, width = 10)


# table studies/time series per area
area <- bio_short %>% 
  group_by(AREA_SQ_KM) %>% 
  summarise(studies = length(unique(STUDY_ID)),
            time_series = length(unique(STUDY_ID_PLOT)))

# table areas over 1km ----
area_over1 <- bio_short %>% 
  group_by(AREA_SQ_KM) %>% 
  filter(AREA_SQ_KM > 1.000000e+00) %>% 
  summarise(studies = length(unique(STUDY_ID)),
            time_series = length(unique(STUDY_ID_PLOT)))
colSums(area_over1)


#area_under96 <- bio_short %>% 
  #group_by(AREA_SQ_KM) %>% 
  #filter(AREA_SQ_KM <= 9.663440e+01) %>% 
  #summarise(studies = length(unique(STUDY_ID)),
            #time_series = length(unique(STUDY_ID_PLOT)))
#colSums(area_under96)


#(area_hist5 <- ggplot(area_under96, aes(x = AREA_SQ_KM)) +
    #geom_histogram()+
    #theme_clean())

# number of studies with permanent plots ----
has_plot <- bio_short %>% 
  filter(AREA_SQ_KM > 1.000000e+00) %>% 
  group_by(HAS_PLOT) %>% 
  summarise(studies = length(unique(STUDY_ID)),
            time_series = length(unique(STUDY_ID_PLOT)))
colSums(has_plot)

# number lat long ----
number_ll <- bio_short %>% 
  filter(AREA_SQ_KM > 1.000000e+00) %>% # studies above 1km2
  filter(HAS_PLOT == "Y") %>% # permanent plots
  group_by(NUMBER_LAT_LONG) %>% 
  summarise(studies = length(unique(STUDY_ID)),
            time_series = length(unique(STUDY_ID_PLOT)))
colSums(number_ll)
# last row does not make sense: more lat/long than time-series

# cent_lat = latitude ----
no_latitude <- bio_short %>% 
  filter(AREA_SQ_KM > 1.000000e+00) %>% 
  #filter(HAS_PLOT == "Y") %>% 
  mutate(test = CENT_LAT == LATITUDE) %>% 
  group_by(test) %>% 
  summarise(studies = length(unique(STUDY_ID)),
            time_series = length(unique(STUDY_ID_PLOT)))

no <- bio_short %>% 
  mutate(test = CENT_LAT == LATITUDE) %>% 
  filter(test == TRUE) %>% 
  filter(AREA_SQ_KM > 1.000000e+00)

# jaccard ~ area ----
# checking distribution
hist(data1$Jtu)
hist(data1$AREA_SQ_KM)

# model
j_a <- lm(Jtu ~ AREA_SQ_KM, data = data1)
summary(j_a)
plot(j_a)

# visualisation
(j_a <- ggplot(data1, aes(x = log(AREA_SQ_KM), y = Jtu)) +
    geom_point(colour = "#483D8B", alpha = 0.3, size = 2) +
    geom_smooth(method = glm, colour = "#483D8B", fill = "#483D8B", alpha = 0.3, size = 2) +
    theme_classic() +
    labs(x = "\nlog Area (km²) ", y = "Jaccard\n"))
         
ggsave(filename = "outputs/jacc_area.png", device = "png", width = 8, height = 6)

# bayesian
j_a_b <- MCMCglmm(Jtu ~ log(AREA_SQ_KM),  data = data1)
summary(j_a_b)



# quantifications ----
acc_75 <- data1 %>% 
  filter(scaleacc_25 > 0.9)
# 96.7% above 0.75
# 90.5% above 0.9

hpd_10 <- data1 %>% 
  filter(scalehpd_25 < 0.25)
# 96.4 below 0.1
# 99.2% below 0.25

taxa_check <- data1 %>% 
  filter(TAXA == "Terrestrial plants")
