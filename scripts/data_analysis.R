# Data transformations
# Dani Gargya, daniela@gargya.de
# March 2020

####### workflow: ----
# calculate jaccard dissimilarity (with/without rarefaction?)
  # filter max and min year only
  # summarise all abundances
  # create site species matrix per plot
  # calculate jaccard dissimilarity
# extract accessibility score/ HPD (what scale?)
# run model
  #create broad grid cell random effect


# calculating jaccard turnover ----
# Load packages ----
library(tidyverse) # (Contains loads of useful functions)
library(ggplot2) # (Package for making nice graphs)
library(vegan)
library(rgdal) # to read in and save spatial data
library(raster) # to allow creation, reading, manip of raster data

# load data ----
# avatar_data <- read.csv("data/avatar_fauna.csv")
# bio <- read.csv("data/bio.csv")
hpd <- raster("data/gpw_v4_population_density_rev11_2015_30_sec.tif")
wp <- raster("data/ppp_2015_1km_Aggregated.tif")


# data manipulation ----
bio_turnover <- bio %>% 
  dplyr::select(STUDY_ID_PLOT, YEAR, GENUS_SPECIES, sum.allrawdata.ABUNDANCE) %>% 
  #group_by(STUDY_ID_PLOT) %>% 
  #filter(YEAR %in% c(max(YEAR), min(YEAR))) %>% 
  #mutate(number_plots = length(unique(YEAR))) %>% 
  #filter(number_plots == 2) %>% 
  filter(STUDY_ID_PLOT %in% c("10_1" , "10_2")) %>% 
  group_by(STUDY_ID_PLOT, YEAR, GENUS_SPECIES) %>% 
  summarise(Abundance = sum(sum.allrawdata.ABUNDANCE)) %>% 
  ungroup()

# spreading into matrix
bio_t_matrix <- bio_turnover %>% 
  #group_by(STUDY_ID_PLOT) %>% 
  spread(GENUS_SPECIES, Abundance, fill = 0) %>% 
  dplyr::select(- STUDY_ID_PLOT, -YEAR)

# calculating jaccard
jaccard1 <- vegdist(bio_t_matrix, method = "jaccard")

# function to calculate jaccard ----
# select each unique STUDY_ID_PLOT
# spread data into matrix
# calculate jaccard
# bring number back to right ID

# create different dataframes
bio_t_list <- split(bio_turnover, unique(bio_turnover$STUDY_ID_PLOT))
#str(bio_t_list)
#bio_t_list2 <- lapply(unique(bio_turnover$STUDY_ID_PLOT), function(x) bio_turnover[bio_turnover$STUDY_ID_PLOT == x,])



# function spreading
#spread.matrix <- function(x, y){
  #tidyr::spread(x,y, fill = 0)
#}

#spread.matrix(x = bio_turnover$GENUS_SPECIES, y = bio_turnover$Abundance)
#matrix <- spread.matrix(x = bio_t_list[[i]]$GENUS_SPECIES, y = bio_t_list[[i]]$Abundance)


# create empty list
jaccard_list <- list()
IDs <- unique(bio_turnover$STUDY_ID_PLOT)

# for loop
#for (i in 1:length(bio_t_list)) {
  #matrix <- tidyr::spread(bio_t_list[[i]]$GENUS_SPECIES, bio_t_list[[i]]$Abundance, fill = 0)
  #matrix_s <- select(matrix, - STUDY_ID_PLOT, -YEAR)
  #jaccard <- vegdist(matrix_s, method = "jaccard")
  #dat <- data.frame(STUDY_ID_PLOT, jaccard)
  #jaccard_list[[i]] <- dat
#}

for (i in 1:length(bio_t_list)) {
  bio_t_list[[i]] <- bio_t_list[[i]] %>% 
    spread(GENUS_SPECIES, Abundance, fill = 0) %>% 
    dplyr::select(-STUDY_ID_PLOT, -YEAR) %>% 
    vegdist(method = "jaccard") -> jaccard
    name <- paste(i, "study_ID",  sep = "_" )
    jaccard_list[[name]] <- jaccard
  }

view(jaccard_list)
# data transformation hpd ----
# explore data
hpd
plot(hpd)
hpd_center <- crop(hpd, extent(-100, 100, -90, 90))
values <- getValues(hpd, 100, 100)
head(values)

plot(hpd_center,
     breaks = c(1, 5, 25, 250, 1000), 
     #col=colorRampPalette(c("blue", "limegreen", "yellow", "darkorange", "red"))(50), # using same colours as in the GEE for the legend
     #legend.args = list(text = 'Annual mean LST (°C)', side = 4, # adding a legend title
                        #font = 1, line = 2.5, cex = 0.75),
     axes=FALSE, # removing the axis
     main = 'Mean land surface temperature in the UK 2017')


# extract values at specific lat/longs
SP <- bio_short %>% 
  dplyr::select(LATITUDE, LONGITUDE)

SP <- SP %>% 
  distinct(LATITUDE, .keep_all = TRUE)

e <- extract(hpd, SP)
e

# scale between 0 and 1
hpd_df <- as.data.frame(hpd)
#hpd_scale <- hpd %>%
  #mutate(scalepop=(pop-min(pop))/(max(pop)-min(pop)))


# wp dataset ----
# remove NA
plot(wp >= 0, wp<= 5000)
plot(wp, col=colorRampPalette(c("blue", "limegreen", "yellow", "darkorange", "red"))(5),
     breaks = c(1, 5, 25, 250, 1000))
head(wp)
minValue(wp)
maxValue(wp)
wp <- setMinMax(wp)
plot(wp)
hist(wp)

# accessibility ----
acc1 <- raster("data/A-0000000000-0000000000.tif")
acc2 <- raster("data/A-0000000000-0000032768.tif")

plot(acc1)
plot(acc2)
acc1
acc2
