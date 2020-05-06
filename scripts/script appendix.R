# Full R script analysing biodiversity changes across levels of accessibility
# Dissertation project 2020, University of Edinburgh
# Daniela Gargya
# May 2020

# Load libraries ----
library(tidyverse)
library(ggplot2) 
library(vegan)
library(betapart)
library(rgdal)
library(raster) 
library(labdsv)
library(dggridR)
library(brms)
library(ggeffects)
library(tidybayes)
library(bayesplot)
library(maps)
library(ggthemes)
library(treemap)
library(RColorBrewer)
library(treemapify)
library(ggpubr)
library(mapdata)
library(gridExtra)

# Load data ----
# BioTIME
biotime_full <- read.csv("data/BioTIMEQuery02_04_2018.csv") 
biotime_meta <- read.csv("data/BioTIMEMetadata_02_04_2018.csv")
# accessed through website: http://biotime.st-andrews.ac.uk/downloadArea.php

# accessibility data
aa <- raster("data/2015_accessibility_to_cities_v1.0/2015_accessibility_to_cities_v1.0.tif")
# accessed through website: https://malariaatlas.org/explorer/#/

# human population density data
hpd <- raster("data/gpw_v4_population_density_rev11_2015_30_sec.tif")
# accessed through website: https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11

### filter biotime data ----
### data inclusion criteria
# terrestrial realm
# minimum study duration 5 years
# at least 15 studies per taxa
# no more than 5000 plots per study
# at least 2 survey points in each plot
# filter by fixed plot

bio <- biotime_meta %>% 
  filter(REALM == "Terrestrial") %>% # only terrestial species
  group_by(STUDY_ID) %>% 
  mutate(duration = END_YEAR - START_YEAR) %>% 
  filter(!duration < 5) %>% # minimum duration of 5 years
  #filter(!DATA_POINTS < 5) %>% 
  group_by(TAXA) %>% 
  mutate(studies_taxa=length(unique(STUDY_ID))) %>% 
  filter(!studies_taxa < 15) %>%  # minimum 15 studies per taxa
  ungroup() %>% 
  # merge datasets
  left_join(biotime_full, by = "STUDY_ID") %>% 
  group_by(STUDY_ID) %>% 
  unite(STUDY_ID_PLOT, STUDY_ID, PLOT, sep = "_", remove=F) %>% 
  group_by(STUDY_ID) %>% 
  mutate(maxplots = length(unique(STUDY_ID_PLOT))) %>% 
  filter(!maxplots > 5000)  %>% 
  ungroup() %>% 
  group_by(STUDY_ID_PLOT) %>% 
  filter(HAS_PLOT == "Y") %>% 
  filter(YEAR %in% c(max(YEAR), min(YEAR))) %>% 
  mutate(number_plots = length(unique(YEAR))) %>% 
  filter(number_plots == 2) %>% # have min and max year per plot only
  dplyr::select(STUDY_ID, PLOT, STUDY_ID_PLOT, START_YEAR, END_YEAR, duration, TAXA, LATITUDE, LONGITUDE, YEAR, sum.allrawdata.ABUNDANCE, GENUS_SPECIES, LATITUDE, LONGITUDE, AREA_SQ_KM, NUMBER_OF_SAMPLES, ABUNDANCE_TYPE, AB_BIO, BIOMASS_TYPE, PROTECTED_AREA) %>% 
  ungroup()

### calculating temporal turnover as component of beta diversity ----

# data manipulation ----
# change biomass and density data to 1 -> absence/presence
bio$sum.allrawdata.ABUNDANCE[bio$sum.allrawdata.ABUNDANCE < 1] <- 1

# summarise abundances per year per species per plot
bio_turnover <- bio %>% 
  dplyr::select(STUDY_ID_PLOT, YEAR, GENUS_SPECIES, sum.allrawdata.ABUNDANCE) %>% 
  group_by(STUDY_ID_PLOT, YEAR, GENUS_SPECIES) %>% 
  summarise(Abundance = sum(sum.allrawdata.ABUNDANCE)) %>% 
  ungroup()

# creating empty dataframe ----
beta_Jaccard <- data.frame(matrix(ncol = 6, nrow = length(unique(bio_turnover$STUDY_ID_PLOT)))) 
names(beta_Jaccard) <- c("STUDY_ID_PLOT", "duration_plot", "richness_change", "Jbeta", "Jtu", "Jne") 
i = 1

# for loop with betapart (code help received from Isla) ----
for (i in 1:length(unique(bio_turnover$STUDY_ID_PLOT))) {
  StudyIDPlot <- as.character(unique(bio_turnover$STUDY_ID_PLOT)[i])
  sub_bio_abundance <- filter(bio_turnover, 
                              STUDY_ID_PLOT == StudyIDPlot)
  duration_plot <- (max(sub_bio_abundance$YEAR) - min(sub_bio_abundance$YEAR)) + 1
  sub_bio_abundance_wider <- pivot_wider(sub_bio_abundance, names_from = GENUS_SPECIES, 
                                         values_from = Abundance, 
                                         values_fill = list(Abundance = 0))
  sub_bio_abundance_matrix <- dplyr::select(sub_bio_abundance_wider, -STUDY_ID_PLOT, -YEAR) 
  sub_bio_presence_matrix <- with(sub_bio_abundance_matrix, ifelse(sub_bio_abundance_matrix > 0,1,0))
  J_components <- beta.pair(sub_bio_presence_matrix, index.family='jaccard')	# distance
  richness_change <- rowSums(sub_bio_presence_matrix)[2] - rowSums(sub_bio_presence_matrix)[1]
  Jbeta <- J_components$beta.jac
  Jtu <- J_components$beta.jtu
  Jne <- J_components$beta.jne
  beta_Jaccard[i,] <- c(StudyIDPlot, duration_plot, richness_change, Jbeta, Jtu, Jne)
  
  i = i+1
}

### extracting values accessibility ----
# df of all unique lat/long values
bio_short <- bio %>% 
  distinct(STUDY_ID_PLOT, .keep_all = TRUE)

SP <- bio_short %>% 
  dplyr::select(LATITUDE, LONGITUDE, STUDY_ID_PLOT) %>% 
  distinct(LATITUDE, .keep_all = TRUE)

# adjusting to crs
points <- cbind(SP$LONGITUDE, SP$LATITUDE)
sppoints <- SpatialPoints(points, proj4string=CRS('+proj=longlat +datum=WGS84'))
tp <- spTransform(sppoints, crs(aa))

# extract numbers at different scales
e <- extract(aa, tp)
e_25 <- extract(aa, tp, buffer = 25000, fun = mean)
e_50 <- extract(aa, tp, buffer = 50000, fun = mean)
e_100 <- extract(aa, tp, buffer = 100000, fun = mean)

bio_aa_2 <- cbind(SP, e, e_2, e_5, e_25, e_50, e_75, e_100)

# drop NA according to scale I am looking at
# scale extracted values between zero and one
bio_aa_short <- bio_aa_2 %>% 
  drop_na(e_25) %>% 
  mutate(scaleacc_25= 1 - ((e_25 -min(e_25))/(max(e_25)-min(e_25))))

# join to full dataset
bio_full_acc <- bio_aa_short %>% 
  right_join(bio_short, by = "LATITUDE") %>% 
  dplyr::select(-STUDY_ID_PLOT.x, -LONGITUDE.y) %>% 
  rename(STUDY_ID_PLOT = STUDY_ID_PLOT.y,
         LONGITUDE = LONGITUDE.x) %>% 
  dplyr::select(STUDY_ID_PLOT, scaleacc_25, e, e_25, e_50, e_100)


### extracting values hpd ----

# turn lat/long values into right CRS format
tp_hpd <- spTransform(sppoints, crs(hpd))

# extract long/lat from raster at different scales
e_hpd <- extract(hpd, tp_hpd)
e_hpd25 <- extract(hpd, tp_hpd, buffer = 25000, fun = mean)
e_hpd50 <- extract(hpd, tp_hpd, buffer = 50000, fun = mean)
e_hpd100 <- extract(hpd, tp_hpd, buffer = 100000, fun = mean)

# bind extracted values to dataframe
bio_hpd2 <- cbind(SP, e_hpd, e_hpd2, e_hpd5, e_hpd25, e_hpd50, e_hpd75, e_hpd100)

# drop NA according to scale I am looking at
# scale extraccted values between zero and one
bio_hpd_short <- bio_hpd2 %>% 
  drop_na(e_hpd25) %>% 
  mutate(scalehpd_25= (e_hpd25-min(e_hpd25))/(max(e_hpd25)-min(e_hpd25)))

# join to full dataset
bio_full_hpd <- bio_hpd_short %>% 
  right_join(bio_short, by = "LATITUDE") %>% 
  dplyr::select(-STUDY_ID_PLOT.x, -LONGITUDE.y) %>% 
  rename(STUDY_ID_PLOT = STUDY_ID_PLOT.y,
         LONGITUDE = LONGITUDE.x) %>% 
  dplyr::select(STUDY_ID_PLOT, scalehpd_1, e_hpd, e_hpd2, e_hpd5, e_hpd25, e_hpd50, e_hpd75, e_hpd100)


### creating global grid cell variable ----
#Construct a global grid with cells approximately 100km² (res= 12 equivalates	95.97785km² cell area)
dggs <- dgconstruct(res = 12, metric=FALSE, resround='down')

#Get the corresponding grid cells for each lat-long pair
bio_short$cell <- dgGEO_to_SEQNUM(dggs,bio_short$LONGITUDE,bio_short$LATITUDE)$seqnum

#Get the number of time-series in each equally-sized cell
biocounts   <- bio_short %>% group_by(cell) %>% summarise(count=n())


### joining dataset with all important variables ----
data1 <- beta_Jaccard %>% 
  left_join(bio_full_acc, by = "STUDY_ID_PLOT") %>% 
  left_join(bio_full_hpd, by = "STUDY_ID_PLOT") %>% 
  left_join(bio_short, by = "STUDY_ID_PLOT")

### Data analysis ----
# Models ----
# main model
mo_tu_simp6 <- brm(bf(Jtu ~ scaleacc_25 + scalehpd_25 + duration_plot + TAXA +
                        (1|cell) + (1|STUDY_ID)),
                   family = zero_one_inflated_beta(), 
                   data = data1,
                   iter = 4000,
                   warmup = 1000,
                   inits = '0',
                   control = list(adapt_delta = 0.85),
                   cores = 4, chains = 4)

save(mo_tu_simp6, file = "outputs/mo_tu_simp6.RData")

# sensitivity analysis plants
data1_plants <- data1 %>% 
  filter(TAXA == "Terrestrial plants")

mo_tu_simp_plants <- brm(bf(Jtu ~ scaleacc_25 + scalehpd_25 + duration_plot +
                              (1|STUDY_ID)), # or with (1|cell)?
                         family = zero_one_inflated_beta(), 
                         data = data1_plants,
                         iter = 4000,
                         warmup = 1000,
                         inits = '0',
                         control = list(adapt_delta = 0.85),
                         cores = 4, chains = 4)

save(mo_tu_simp_plants, file = "outputs/mo_tu_simp_plants.RData")

# sensitivity analysis scale (1km)
mo_tu_simp8 <- brm(bf(Jtu ~ scaleacc_1 + scalehpd_1 + duration_plot + TAXA +
                        (1|cell) + (1|STUDY_ID)), 
                   family = zero_one_inflated_beta(), 
                   data = data1,
                   iter = 4000,
                   warmup = 1000,
                   inits = '0',
                   control = list(adapt_delta = 0.85),
                   cores = 4, chains = 4)

save(mo_tu_simp8, file = "outputs/mo_tu_simp8.RData")


# sensitivity analysis scale (50km)
mo_tu_simp7 <- brm(bf(Jtu ~ scaleacc_50 + scalehpd_50 + duration_plot + TAXA +
                        (1|STUDY_ID)), # or with (1|cell)?
                   family = zero_one_inflated_beta(), 
                   data = data1,
                   iter = 4000,
                   warmup = 1000,
                   inits = '0',
                   control = list(adapt_delta = 0.85),
                   cores = 4, chains = 4)

save(mo_tu_simp7, file = "outputs/mo_tu_simp7.RData")

# sensitivity analysis scale 100km
mo_tu_scale100 <- brm(bf(Jtu ~ scaleacc_100 + scalehpd_100 + duration_plot + TAXA +
                           (1|cell) + (1|STUDY_ID)),
                      family = zero_one_inflated_beta(), 
                      data = data1,
                      iter = 4000,
                      warmup = 1000,
                      inits = '0',
                      control = list(adapt_delta = 0.85),
                      cores = 2, chains = 4)

save(mo_tu_scale100, file = "outputs/mo_tu_scale100.RData")

#### Data visualisation ----
# create colour palette ----
display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n = 8, name = "Dark2")
display.brewer.pal(n=4, name = "Set1")
brewer.pal(n=4, name = "Set1")
taxa.palette <- c("#D95F02", "#7570B3", "#E7298A", "#E6AB02")
names(taxa.palette) <- levels(data1$TAXA)

# clean theme ----
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
          legend.title = element_text(size = 12, face = "bold"),                              
          legend.position = c(0.2, 0.8))
}



#### panel time-series ----
# spatial distribution of time-series ----
(map_studies2 <- ggplot(bio_short,
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
taxa_sum <- bio_short %>%  group_by(TAXA) %>% tally

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

#### RQ1: turnover ~ accessibility + duration ----
# prediction
pred_acc_d <- ggpredict(mo_tu_simp6, terms = c("scaleacc_25", "duration_plot"))

pred_acc_d$Monitoring <- factor(pred_acc_d$group, levels = c("6", "19.1", "32.2"),
                                labels = c("Short-term", "Moderate", "Long-term"))

# calculate change
pred_acc <- ggpredict(mo_tu_simp6, terms = c("scaleacc_25"))

# model vis
(graph_acc <- ggplot() +
    geom_line(data = pred_acc_d, aes(x = x, y = predicted, color = Monitoring),
              size = 2) +
    geom_ribbon(data = pred_acc_d, aes(ymin = conf.low, ymax = conf.high, 
                                       x = x, fill = Monitoring), alpha = 0.1) +
    geom_point(data = data1, aes(x = scaleacc_25, y = Jtu),
               alpha = 0.1, size = 2) +
    theme_clean() +
    scale_fill_manual(values = c("#7CFC00", "#66A61E", "#1B9E77")) +
    scale_colour_manual(values = c("#7CFC00", "#66A61E", "#1B9E77")) +
    labs(x = "\nAccessibility (proportion)", y = "Turnover\n") +
    theme(legend.position = c(0.2, 0.8)))

ggsave(graph_acc, file = "outputs/graph_acc.png", width = 7, height = 5)

(acc_d <- ggMarginal(graph_acc, type="density", size = 3, fill = "#666666", xparams = list(fill = "#A6761D")))

ggsave(acc_d, file = "outputs/graph_acc_den.png" , width = 7, height = 5)

#### RQ2: turnover ~ hpd ----
# prediction
pred_hpd <- ggpredict(mo_tu_simp6, terms = c("scalehpd_25"))

# model vis
(graph_hpd <- ggplot() +
    geom_line(data = pred_hpd, aes(x = x, y = predicted),
              size = 2, color = "#A6761D") +
    geom_ribbon(data = pred_hpd, aes(ymin = conf.low, ymax = conf.high, 
                                     x = x), fill = "#A6761D", alpha = 0.1) +
    geom_point(data = data1, aes(x = scalehpd_25, y = Jtu),
               alpha = 0.1, size = 2) +
    theme_clean() +
    labs(x = "\nHuman population density (proportion)", y = "Turnover\n") +
    theme(legend.position = c(0.85, 0.8)))

ggsave(graph_hpd, file = "outputs/graph_hpd.png", width = 7, height = 5)

(hpd_d <- ggMarginal(graph_hpd, type="density", size = 3, fill = "#666666", xparams = list(fill = "#A6761D")))

ggsave(hpd_d, file = "outputs/graph_hpd_den.png" , width = 7, height = 5)

#### RQ3: turnover ~ taxa ----
# raincloud plot
# by Ben Marwick
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

# raincloud plot taxa jtu
(raincloud_taxa <- 
    ggplot(data = data1, 
           aes(x = reorder(TAXA, desc(Jtu)), y = Jtu, fill = TAXA)) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.8) +
    geom_point(aes(y = Jtu, color = TAXA), 
               position = position_jitter(width = 0.15), size = 1, alpha = 0.1) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
    labs(y = "\n Turnover", x = NULL) +
    guides(fill = FALSE, color = FALSE) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = taxa.palette) +
    scale_colour_manual(values = taxa.palette) +
    coord_flip() +
    theme_clean())

ggsave(raincloud_taxa, filename = "outputs/raincloud_taxa_graph.png",
       height = 5, width = 8)

# effect size graph
pred_taxa <- ggpredict(mo_tu_simp6, terms = c("TAXA"))

(graph_taxa <- ggplot(pred_taxa, aes(x = reorder(x, desc(predicted)), y = predicted, color = x)) +
    geom_point(aes(size = 1)) +
    geom_pointrange(aes(ymin = conf.low, ymax =  conf.high), size = 2) +
    scale_fill_manual(values = taxa.palette) +
    scale_colour_manual(values = taxa.palette) +
    theme_clean() +
    theme(legend.position = "none",
          axis.text.y=element_blank(),
          axis.title.y=element_blank()) +
    geom_hline(yintercept=0, linetype="dashed", size=1) +
    labs(x = NULL, y = "Effect size\n") +
    coord_flip())

ggsave(graph_taxa, filename = "outputs/taxa_model.png",
       height = 5, width = 5)

panel_taxa <- ggarrange(raincloud_taxa, graph_taxa, labels = c("A", "B"), ncol = 2)
panel_taxa2 <- grid.arrange(raincloud_taxa, graph_taxa, respect=TRUE, ncol = 2,  widths=c(1.5,1))

#graph_taxa = graph_taxa + theme(aspect.ratio=2)

#p <- grid.arrange(raincloud_taxa, graph_taxa,  ncol = 2) # one is square, the other thinner  


ggsave(p, filename ="outputs/panel_taxa3.png",
       height = 5, width = 8)

#### sensitivtiy analysis ----
# plants vs all ----
pred_plants <- ggpredict(mo_tu_simp_plants, terms = c("scaleacc_25"))
pred_all <- ggpredict(mo_tu_simp6, terms = c("scaleacc_25"))

data1_plants <- data1 %>% 
  filter(TAXA == "Terrestrial plants")

data1_nopl <- data1 %>% 
  filter(!TAXA == "Terrestrial plants")

(graph_plants <- ggplot() +
    geom_line(data = pred_plants, color = "#E6AB02", aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_plants, color = "#E6AB02", aes(ymin = conf.low, ymax = conf.high, 
                                                           x = x), alpha = 0.1, fill = "#E6AB02") +
    geom_point(data = data1_plants, color = "#E6AB02", aes(x = scaleacc_25, y = Jtu),
               alpha = 0.1, size = 2) +
    geom_line(data = pred_all, color = "#66A61E", aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_all, color = "#66A61E", aes(ymin = conf.low, ymax = conf.high, 
                                                        x = x), alpha = 0.1, fill = "#66A61E") +
    geom_point(data = data1_nopl, color = "#66A61E", aes(x = scaleacc_25, y = Jtu),
               alpha = 0.1, size = 2) +
    theme_clean() +
    labs(x = "\nAccessibility", y = "Turnover\n"))

ggsave(graph_plants, filename = "outputs/graph_plants.png",  height = 5, width = 8)

# scales ----
pred_1 <- ggpredict(mo_tu_simp8, terms = c("scaleacc_1"))
pred_25 <- ggpredict(mo_tu_simp6, terms = c("scaleacc_25"))
pred_50 <- ggpredict(mo_tu_simp7, terms = c("scaleacc_50"))
pred_100 <- ggpredict(mo_tu_scale100, terms = c("scaleacc_100"))

(graph_scale <- ggplot() +
    geom_line(data = pred_1, color = "#E41A1C", aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_1, color = "#E41A1C", aes(ymin = conf.low, ymax = conf.high, 
                                                      x = x), alpha = 0.1, fill = "#E41A1C") +
    geom_line(data = pred_25, color = "#377EB8" , aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_25, color = "#377EB8" , aes(ymin = conf.low, ymax = conf.high, 
                                                        x = x), alpha = 0.1, fill = "#377EB8" ) +
    geom_line(data = pred_50, color = "#4DAF4A" , aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_50, color = "#4DAF4A", aes(ymin = conf.low, ymax = conf.high, 
                                                       x = x), alpha = 0.1, fill = "#4DAF4A" ) +
    geom_line(data = pred_100, color = "#984EA3", aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_100, color = "#984EA3", aes(ymin = conf.low, ymax = conf.high, 
                                                        x = x), alpha = 0.1, fill = "#984EA3") +
    theme_clean() +
    labs(x = "\nAccessibility", y = "Turnover\n") +
    ylim(0, 1))

ggsave(graph_scale, filename = "outputs/graph_scale.png", height = 5, width = 8)
#### other additional graphs ----
# raincloud plot taxa accessibility ----
(raincloud_acc <- 
   ggplot(data = data1, 
          aes(x = reorder(TAXA, desc(scaleacc_25)), y = scaleacc_25, fill = TAXA)) +
   geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.8) +
   geom_point(aes(y = scaleacc_25, color = TAXA), 
              position = position_jitter(width = 0.15), size = 1, alpha = 0.1) +
   geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
   labs(y = "\nAccessibility (proportion)", x = NULL) +
   guides(fill = FALSE, color = FALSE) +
   scale_fill_manual(values = taxa.palette) +
   scale_colour_manual(values = taxa.palette) +
   coord_flip() +
   theme_clean())

ggsave(raincloud_acc, filename = "outputs/raincloud_taxa_acc_graph.png",
       height = 5, width = 8)

# protected areas ----
(graph_proteccted_areas <- ggplot() +
   geom_point(data = data1, aes(x = scaleacc_25, y = Jtu, colour = PROTECTED_AREA),
              alpha = 0.8, size = 1.5) +
   scale_colour_manual(values = c("#E6AB02", "#666666")) +
   theme_clean() +
   theme(legend.position = "none") +
   labs(x = "\nAccessibility", y = "Turnover\n"))

ggsave(graph_proteccted_areas, filename = "outputs/graph_protected_areas.png",
       height = 5, width = 8)

# graph duration and turnover ----
pred_duration <- ggpredict(mo_tu_simp6, terms = c("duration_plot"))

(graph_duration <- ggplot() +
    geom_line(data = pred_duration, aes(x = x, y = predicted),
              size = 2, color = "#666666") +
    geom_ribbon(data = pred_duration, aes(ymin = conf.low, ymax = conf.high, 
                                          x = x), fill = "#666666", alpha = 0.1) +
    theme_clean() +
    labs(x = "\nMonitoring duration (in years)", y = "Turnover\n") +
    ylim(0, 1))

ggsave(graph_duration, file = "outputs/graph_duration.png", width = 7, height = 5)

# histogram distribution of random and real values ----
# accessibility
bio_short <- bio %>% 
  distinct(STUDY_ID_PLOT, .keep_all = TRUE)

SP <- bio_short %>% 
  dplyr::select(LATITUDE, LONGITUDE, STUDY_ID_PLOT) %>% 
  distinct(LATITUDE, .keep_all = TRUE)

library(generator)
fake_lat <- r_latitudes(1023)
fake_long <- r_longitudes(1023)

fake_ll <- SP %>% 
  mutate(fake_lat = c(fake_lat),
         fake_long = c(fake_long)) %>% 
  dplyr::select(- LATITUDE, -LONGITUDE)

f_points <- cbind(fake_ll$fake_long, fake_ll$fake_lat)


f_sppoints <- SpatialPoints(f_points, proj4string=CRS('+proj=longlat +datum=WGS84'))
f_tp <- spTransform(f_sppoints, crs(aa))

f_e <- extract(aa, f_tp)

f_bio_aa <- cbind(fake_ll, f_e)
f_bio_aa_short <- na.omit(f_bio_aa)

# hpd
f_tp_hpd <- spTransform(f_sppoints, crs(hpd))

f_e_hpd <- extract(hpd, f_tp_hpd)

f_bio_hpd <- cbind(fake_ll, f_e_hpd)
f_bio_hpd_short <- na.omit(f_bio_hpd)

# figure overall
png("outputs/fake_ll.png", width=600, height=350)
par(
  mfrow=c(1,2),
  mar=c(4,4,1,0)
)
hist(f_bio_aa_short$f_e, col=rgb(1,0,0,0.5), xlab="Accessibility score", ylab ="Number of points", main = "") # more normal distributed?
hist(data1$e, col=rgb(0,0,1,0.5), add=T)
hist(f_bio_hpd_short$f_e_hpd, col=rgb(1,0,0,0.5), xlab="Human population density score", ylab ="Number of points", main = "")
hist(data1$e_hpd, col=rgb(0,0,1,0.5), add=T)
legend("topright", legend=c("Random","BioTIME"), col=c(rgb(1,0,0,0.5), 
                                                       rgb(0,0,1,0.5)), pt.cex=2, pch=15)
dev.off()

# histogram acc and hpd ----
png("outputs/hist_acc_hpd.png", width=600, height=350)
hist(data1$scaleacc_25, col = "#1B9E77", xlab="Accessibility/Human population density", ylab ="Number of points", main = "")
hist(data1$scalehpd_25, col = "#A6761D", add =T)
legend("top", legend=c("Accessibility","Human population density"), col=c("#1B9E77", "#A6761D"), pt.cex=2, pch=15)
dev.off()