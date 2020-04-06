# Data visualisation
# Dani Gargya
# March 2020

# Load data ----
biotime_full <- read.csv("data/BioTIMEQuery02_04_2018.csv")
biotime_meta <- read.csv("data/BioTIMEMetadata_02_04_2018.csv")
bio <- read.csv("data/bio.csv")

bio_short <- bio %>% 
  distinct(STUDY_ID_PLOT, STUDY_ID, CENT_LAT, CENT_LONG, TAXA, TOTAL, START_YEAR, END_YEAR, NUMBER_OF_SAMPLES, AREA_SQ_KM, HAS_PLOT, HAS_PLOT, NUMBER_LAT_LONG)

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
library(dggridR)
library(devtools)

#find_rtools()
#pkgbuild::has_rtools

#install_url('https://cran.r-project.org/src/contrib/Archive/dggridR/dggridR_2.0.3.tar.gz')

packageurl <- 'https://cran.r-project.org/src/contrib/Archive/dggridR/dggridR_2.0.3.tar.gz'
install.packages(packageurl, repos=NULL, type="source")

install_github('r-barnes/dggridR', vignette=TRUE)
install_github('r-barnes/dggridR', vignette=TRUE)

devtools::install_github("karthik/wesanderson")


# RQ1: jaccard ~ accessibility ----
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


ggplot() +
  geom_line(data = predictions, aes(x = x, y = predicted),
            size = 2) +
  geom_ribbon(data = predictions, aes(ymin = conf.low, ymax = conf.high, 
                                      x = x), alpha = 0.1) +
  geom_point(data = data1, aes(x = scaleacc, y = Jtu),
             alpha = 0.1, size = 2) +
  #annotate("text", x = -0.65, y = 5, label = "Slope = -0.06, Std. error = 0.01") +  
  #scale_x_continuous(limits = c (0.8, 1)) +
  theme_clean() +
  #scale_fill_manual(values = c("darksalmon", "firebrick3", "firebrick4")) +
  #scale_colour_manual(values = c("darksalmon", "firebrick3", "firebrick4")) +
  labs(x = "\nAccessibility", y = "Jaccard turnover\n")


# RQ2: taxa making raincloud plot ----
# creating theme----
theme_niwot <- function(){
  theme_bw() +
    theme(text = element_text(family = "Helvetica Light"),
          axis.text = element_text(size = 16), 
          axis.title = element_text(size = 18),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 18, vjust = 1, hjust = 0),
          legend.text = element_text(size = 12),          
          legend.title = element_blank(),                              
          legend.position = c(0.95, 0.15), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}

# We will use a function by Ben Marwick
# This code loads the function in the working environment
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

# visualising
(violin_taxa <- ggplot(data1, aes(x = TAXA, y = Jtu)) +
    geom_violin())

# raincloud plot ----
(raincloud_taxa <- 
    ggplot(data = data1, 
           aes(x = reorder(TAXA, desc(Jtu)), y = Jtu, fill = TAXA)) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.8) +
    geom_point(aes(y = Jtu, color = TAXA), 
               position = position_jitter(width = 0.15), size = 1, alpha = 0.1) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
    labs(y = "\nJaccard dissimilarity index", x = NULL) +
    guides(fill = FALSE, color = FALSE) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C")) +
    coord_flip() +
    theme_niwot())

ggsave(raincloud_taxa, filename = "outputs/raincloud_taxa_graph.png",
       height = 5, width = 8)

# RQ3: accessibility and hpd
