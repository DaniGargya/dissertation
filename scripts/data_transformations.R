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


# Load libraries ----
library(tidyverse) # (Contains loads of useful functions)
library(ggplot2) # (Package for making nice graphs)
library(vegan)
library(betapart)
library(rgdal) # to read in and save spatial data
library(raster) # to allow creation, reading, manip of raster data
library(labdsv)

# load data ----
# biodiversity data
# bio <- read.csv("data/bio.csv")

# human population density data
hpd <- raster("data/gpw_v4_population_density_rev11_2015_30_sec.tif")
wp <- raster("data/ppp_2015_1km_Aggregated.tif")

# accessibility data
acc1 <- raster("data/A-0000000000-0000000000.tif")
acc2 <- raster("data/A-0000000000-0000032768.tif")


# calculating jaccard data manipulation ----
bio_turnover <- bio %>% 
  dplyr::select(STUDY_ID_PLOT, YEAR, GENUS_SPECIES, sum.allrawdata.ABUNDANCE) %>% 
  #group_by(STUDY_ID_PLOT) %>% 
  #filter(YEAR %in% c(max(YEAR), min(YEAR))) %>% 
  #mutate(number_plots = length(unique(YEAR))) %>% 
  #filter(number_plots == 2) %>% 
  filter(STUDY_ID_PLOT %in% c("10_1")) %>% 
  group_by(STUDY_ID_PLOT, YEAR, GENUS_SPECIES) %>% 
  summarise(Abundance = sum(sum.allrawdata.ABUNDANCE)) %>% 
  ungroup()

# calculation for 1 ----
# spreading into matrix
bio_t_matrix <- bio_turnover %>% 
  spread(GENUS_SPECIES, Abundance, fill = 0) %>% 
  dplyr::select(- STUDY_ID_PLOT, -YEAR)

# calculating jaccard
jaccard1 <- vegdist(bio_t_matrix, method = "jaccard")

# betapart requires presence/absence matrix for Jaccard calculations of turnover/nestedness
bio_t_matrix_binary <- with(bio_t_matrix, ifelse(bio_t_matrix > 0, 1, 0)) 


# calculation for all ----
# create list for each plot ----
bio_t_list <- split(bio_turnover, unique(bio_turnover$STUDY_ID_PLOT))

# create empty dataframe ----
jaccard_list <- data.frame(STUDY_ID_PLOT = unique(bio_turnover$STUDY_ID_PLOT), jaccard = NA)

# for loop with vegdist ----
for (i in 1:length(bio_t_list)) {
  jaccard_df[[i]] <- bio_t_list[[i]] %>% 
  spread(GENUS_SPECIES, Abundance) %>% 
  dplyr::select(-STUDY_ID_PLOT, -YEAR) %>% 
  vegdist(method = "jaccard", binary = TRUE)
     
  jaccard_list[i, "jaccard"] <- jaccard_df
}

# for loop with betapart ----
for (i in 1:length(bio_t_list)) {
  bio_t_list[[i]] <- bio_t_list[[i]] %>% 
    spread(GENUS_SPECIES, Abundance, fill = 0) %>% 
    dplyr::select(-STUDY_ID_PLOT, -YEAR) -> comm
  
    comm_binary <- with(comm, ifelse(comm > 0, 1, 0)) 
    
    j_components <- beta.pair(comm_binary, index.family = "jaccard")
    
    jtu <- as.matrix (j_components$beta.jtu)
  
  jaccard_list[i, "jaccard"] <- jtu
}
# s change workshop code ----
# calculating between year similarities (NOT DISTANCE!) with Jaccard
Jacsim <- as.matrix(1-vegdist(bio_t_matrix, method='jaccard', binary=TRUE))

J_components <- beta.pair(bio_t_matrix_binary, index.family='jaccard')
Jbeta <- as.matrix(J_components$beta.jac)
Jtu <- as.matrix(J_components$beta.jtu)    

n <- length(unique(bio_turnover$YEAR))

# initialise matrices for calculating turnover
simbaseline <- data.frame(array(NA, dim=c(length(unique(bio_turnover$YEAR)), 4)))
names(simbaseline)<-c('YEAR', 'Jaccard_base', 'Jbeta_base', 'Jtu_base')

simnext <- data.frame(array(NA, dim=c(length(unique(bio_turnover$YEAR)), 4)))
names(simnext)<-c('YEAR', 'Jaccard_base', 'Jbeta_base', 'Jtu_base')

simhind <- data.frame(array(NA, dim=c(length(unique(bio_turnover$YEAR)), 4)))
names(simhind)<-c('YEAR', 'Jaccard_base', 'Jbeta_base', 'Jtu_base')


counter2 <- 1

# baseline
simbaseline[counter2:(counter2+n-2),] <- cbind(
  unique(bio_turnover$YEAR)[2:n],
  Jacsim[2:n],
  Jbeta[2:n],
  Jtu[2:n])

# How consecutive is calculated.
simnext[counter2:(counter2+n-2),] <- cbind(
  unique(bio_turnover$YEAR)[2:n],
  Jacsim[row(Jacsim)-col(Jacsim)==1],
  Jbeta[row(Jbeta)-col(Jbeta)==1],
  Jtu[row(Jtu)-col(Jtu)==1])

# How hindcasting is calculated.  
simhind[counter2:(counter2+n-2),] <- cbind(
  unique(bio_turnover$YEAR)[1:(n-1)],
  Jacsim[row(Jacsim)%in%1:(max(row(Jacsim))-1) & col(Jacsim)==max(col(Jacsim))], 
  Jtu[row(Jtu)%in%1:(max(row(Jtu))-1) & col(Jtu)==max(col(Jtu))])

# combine univariate and turnover metrics
biochange_metrics <- simbaseline[-length(unique(bio_turnover$YEAR)),]
biochange_metrics <- full_join(biochange_metrics, simnext[-length(unique(bio_turnover$YEAR)),], by=c('YEAR'))
biochange_metrics <- full_join(biochange_metrics, simhind[-length(unique(bio_turnover$YEAR)),], by=c('YEAR'))


##	initialise df to store all biochange metrics 
rarefied_metrics <- data.frame()


# add to dataframe for all studies
rarefied_metrics <- bind_rows(rarefied_metrics, biochange_metrics)

rarefied_metrics <- as_tibble(rarefied_metrics) 
# combine with the new metadata
rarefied_metrics <- inner_join(bio_short, rarefied_metrics, by='STUDY_ID_PLOT') 
return(rarefied_metrics)


# data transformation hpd ----
# explore data
hpd
plot(hpd)

hpd_center <- crop(hpd, extent(-100, 100, -90, 90))
values <- getValues(hpd, 100, 100)
head(values)

#hpd_df <- as.data.frame(hpd) # data too big

# create lat/long matrix
SP <- bio_short %>% 
  dplyr::select(LATITUDE, LONGITUDE, STUDY_ID_PLOT) %>% 
  distinct(LATITUDE, .keep_all = TRUE)

# turn lat/long values into right CRS format
points <- cbind(SP$LONGITUDE, SP$LATITUDE)
sppoints <- SpatialPoints(points, proj4string=CRS('+proj=longlat +datum=WGS84'))
tp_hpd <- spTransform(sppoints, crs(hpd))

# extract long/lat from raster
e_hpd <- extract(hpd, tp_hpd)

# bind extracted values to dataframe
bio_hpd <- cbind(SP, e_hpd)

# omit NAs
bio_hpd_short <- na.omit(bio_hpd)

# scale workd population extracted
bio_hpd_scale <- bio_hpd_short %>%
  mutate(scalehpd=(e_hpd-min(e_hpd))/(max(e_hpd)-min(e_hpd)))

# check histogram of values
hist(bio_hpd_scale$scalehpd)
hist(log(bio_hpd_scale$scalehpd))

# wp dataset ----
# data exploration
# remove NAs?

wp
plot(wp)

# plot differently
plot(wp >= 0, wp<= 5000)
plot(wp, col=colorRampPalette(c("blue", "limegreen", "yellow", "darkorange", "red"))(5),
     breaks = c(1, 5, 25, 250, 1000))
head(wp)
minValue(wp)
maxValue(wp)
wp <- setMinMax(wp)
hist(wp)

# create lat/long matrix
SP <- bio_short %>% 
  dplyr::select(LATITUDE, LONGITUDE, STUDY_ID_PLOT) %>% 
  distinct(LATITUDE, .keep_all = TRUE)

# turn lat/long values into right CRS format
points <- cbind(SP$LONGITUDE, SP$LATITUDE)
sppoints <- SpatialPoints(points, proj4string=CRS('+proj=longlat +datum=WGS84'))
tp_wp <- spTransform(sppoints, crs(wp))

# extract long/lat from raster
e_wp <- extract(wp, tp_wp)

# bind extracted values to dataframe
bio_wp <- cbind(SP, e_wp)

# omit NAs
bio_wp_short <- na.omit(bio_wp)

# scale workd population extracted
bio_wp_scale <- bio_wp_short %>%
  mutate(scalewp=(e_wp-min(e_wp))/(max(e_wp)-min(e_wp)))

# check histogram of values
hist(bio_wp_scale$scalewp)
hist(log(bio_wp_scale$scalewp))

# comparing values hpd and wp ----
hpd_wp <- left_join(bio_hpd_scale, bio_wp_scale, by = "STUDY_ID_PLOT")

hpd_wp <- bio_hpd_scale %>% 
  left_join(bio_wp_scale, by = "STUDY_ID_PLOT") %>% 
  dplyr::select(STUDY_ID_PLOT, scalehpd, scalewp)

# accessibility ----
plot(acc1)
plot(acc2)
acc1
acc2

acc <- merge(acc1, acc2)
plot(acc)

acc_df <- as.data.frame(acc)


# extract values at specific lat/longs
SP <- bio_short %>% 
  dplyr::select(LATITUDE, LONGITUDE, STUDY_ID_PLOT) %>% 
  distinct(LATITUDE, .keep_all = TRUE)

points <- cbind(SP$LONGITUDE, SP$LATITUDE)
sppoints <- SpatialPoints(points, proj4string=CRS('+proj=longlat +datum=WGS84'))
tp <- spTransform(sppoints, crs(acc))

e <- extract(acc, tp)

bio_ll <- distinct(bio_short$LATITUDE, .keep_all = TRUE)

bio_acc <- cbind(SP, e)
bio_acc_short <- na.omit(bio_acc)

bio_acc_scale <- bio_acc_short %>%
  mutate(scaleacc=(e-min(e))/(max(e)-min(e)))

hist(log(bio_acc_scale$scaleacc))
