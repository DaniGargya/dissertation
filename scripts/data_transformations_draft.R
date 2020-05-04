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
library(dggridR)

library(sp)
library(rgdal)
library(raster)
library(ggplot2)
library(viridis)
library(rasterVis)

# load data ----
# biodiversity data
# bio <- read.csv("data/bio.csv")

# human population density data
hpd <- raster("data/gpw_v4_population_density_rev11_2015_30_sec.tif")
# https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11
wp <- raster("data/ppp_2015_1km_Aggregated.tif")
# https://www.worldpop.org/geodata/summary?id=24772

# accessibility data (GEE)
acc1 <- raster("data/A-0000000000-0000000000.tif")
acc2 <- raster("data/A-0000000000-0000032768.tif")

aa <- raster("data/2015_accessibility_to_cities_v1.0/2015_accessibility_to_cities_v1.0.tif")
# https://malariaatlas.org/explorer/#/


# calculating jaccard data manipulation ----
bio_turnover3 <- bio %>% 
  dplyr::select(STUDY_ID_PLOT, YEAR, GENUS_SPECIES, sum.allrawdata.ABUNDANCE) %>% 
  #group_by(STUDY_ID_PLOT) %>% 
  #filter(YEAR %in% c(max(YEAR), min(YEAR))) %>% 
  #mutate(number_plots = length(unique(YEAR))) %>% 
  #filter(number_plots == 2) %>% 
  #filter(STUDY_ID_PLOT %in% c("10_1", "10_2", "10_9")) %>% 
  group_by(STUDY_ID_PLOT, YEAR, GENUS_SPECIES) %>% 
  summarise(Abundance = sum(sum.allrawdata.ABUNDANCE)) %>% 
  ungroup()

write.csv(bio_turnover2, "data/bio_turnover2.csv")


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
for (i in seq_along(bio_t_list)) {
  jaccard_df <- bio_t_list[[i]] %>% 
  spread(GENUS_SPECIES, Abundance, fill = 0) %>% 
  dplyr::select(-STUDY_ID_PLOT, -YEAR) %>% 
  vegdist(method = "jaccard", binary = TRUE)
     
  #jaccard_list[i, "jaccard"] <- jaccard_df
}

jacc <- bio_turnover %>% 
  spread(GENUS_SPECIES, Abundance, fill =0) %>% 
  dplyr::select(-YEAR, -STUDY_ID_PLOT)
  
comm_binary <- with(jacc, ifelse(jacc >0,1,0))

j_components <- beta.pair(comm_binary, index.family = "jaccard")

jtu <- as.matrix (j_components$beta.jtu)


# for loop with betapart ----
for (i in seq_along(bio_t_list)) {
  beta_df <- bio_t_list[[i]] %>% 
    spread(GENUS_SPECIES, Abundance, fill = 0)
    dplyr::select(-STUDY_ID_PLOT, -YEAR) -> comm
  
    comm_binary <- with(comm, ifelse(comm > 0, 1, 0)) 
    
    j_components <- beta.pair(comm_binary, index.family = "jaccard")
    
    jtu <- as.matrix (j_components$beta.jtu)
  
  #jaccard_list[i, "jaccard"] <- jtu[1,2]
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
# 1023 unqiue locations

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

# join hpd to all studies
bio_full_hpd <- bio_hpd %>% 
  left_join(bio_hpd_scale, by = "STUDY_ID_PLOT") %>% 
  dplyr::select(-LATITUDE.y, -LONGITUDE.y, -e_hpd.y) %>% 
  rename(LATITUDE = LATITUDE.x) %>% 
  right_join(bio_short, by = "LATITUDE")


# wp dataset ----
# data exploration
# remove NAs?

wp
plot(wp)

# plot differently
plot(wp >= 0, wp<= 5000)
plot(wp, col=colorRampPalette(c("blue", "limegreen", "yellow", "darkorange", "red"))(5),
     breaks = c(1, 5, 25, 250, 1000))

wp[wp > 1200] <- NA
plot(wp)
# mapview?

#image(wp, col= viridis_pal(option="D")(5))
#gplot(wp) +
  #geom_raster(aes(x = x, y = y, fill = value)) +
  # value is the specific value (of reflectance) each pixel is associated with
  #scale_fill_viridis_c() +
  #coord_quickmap() +
  #ggtitle("West of Loch tay, raster plot") +
  #xlab("Longitude") +
  #ylab("Latitude") +
  #theme_classic() +   					    # removes defalut grey background
  #theme(plot.title = element_text(hjust = 0.5),             # centres plot title
   #     text = element_text(size=20),		       	    # font size
    #    axis.text.x = element_text(angle = 90, hjust = 1))  # rotates x axis text

wwp <- setMinMax(wp)
minValue(wwp)
maxValue(wwp)

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

#hist(acc)

# extract values at specific lat/longs
bio_short <- bio %>% 
  distinct(STUDY_ID_PLOT, .keep_all = TRUE)

SP <- bio_short %>% 
  dplyr::select(LATITUDE, LONGITUDE, STUDY_ID_PLOT) %>% 
  distinct(LATITUDE, .keep_all = TRUE)

points <- cbind(SP$LONGITUDE, SP$LATITUDE)
sppoints <- SpatialPoints(points, proj4string=CRS('+proj=longlat +datum=WGS84'))
tp <- spTransform(sppoints, crs(acc))

e <- extract(acc, tp)

bio_acc <- cbind(SP, e)
bio_acc_short <- na.omit(bio_acc)

bio_acc_scale <- bio_acc_short %>%
  mutate(scaleacc=(e-min(e))/(max(e)-min(e)))

hist(log(bio_acc_scale$scaleacc))

# combining hpd and accessibility ----
hpd_acc <- left_join(bio_hpd_scale, bio_acc_scale, by = "STUDY_ID_PLOT")

hpd_acc <- bio_hpd_scale %>% 
  left_join(bio_acc_scale, by = "STUDY_ID_PLOT") %>% 
  dplyr::select(STUDY_ID_PLOT, scalehpd, scaleacc)

# using dataset from website malariaatlas ----
plot(aa)
aa

hist(aa)

points <- cbind(SP$LONGITUDE, SP$LATITUDE)
sppoints <- SpatialPoints(points, proj4string=CRS('+proj=longlat +datum=WGS84'))
tp <- spTransform(sppoints, crs(aa))

e <- extract(aa, tp, buffer = 2000, fun = mean)

bio_aa <- cbind(SP, e)
bio_aa_short <- na.omit(bio_aa)

bio_aa_scale <- bio_aa_short %>%
  mutate(scaleacc=(e-min(e))/(max(e)-min(e)))

hist(log(bio_aa_scale$scaleacc))

bio_full_acc <- bio_aa %>% 
  left_join(bio_aa_scale, by = "STUDY_ID_PLOT") %>% 
  dplyr::select(-LATITUDE.y, -LONGITUDE.y, -e.y) %>% 
  rename(LATITUDE = LATITUDE.x) %>% 
  right_join(bio_short, by = "LATITUDE")

# check how many NAs
bio_full_acc_s <- bio_full_acc %>% 
  dplyr::select(scaleacc, STUDY_ID_PLOT.y) %>% 
  drop_na()

bio_full_hpd_s <- bio_full_hpd %>% 
  drop_na()

# trying to reproject
library(sf)

st_is_longlat(aa)


# create global grid cell ----
#Construct a global grid with cells approximately 1000 miles across
dggs <- dgconstruct(res=12, metric=FALSE, resround='down')

#Get the corresponding grid cells for each earthquake epicenter (lat-long pair)
bio_short$cell <- dgGEO_to_SEQNUM(dggs,bio_short$LONGITUDE, bio_short$LATITUDE)$seqnum

#Converting SEQNUM to GEO gives the center coordinates of the cells
cellcenters   <- dgSEQNUM_to_GEO(dggs,bio_short$cell)

#Get the number of earthquakes in each cell
biocounts   <- bio_short %>% group_by(cell) %>% summarise(count=n())

#Get the grid cell boundaries for cells which had quakes
grid   <- dgcellstogrid(dggs,biocounts$cell,frame=TRUE,wrapcells=TRUE)

#Update the grid cells' properties to include the number of earthquakes
#in each cell
grid <- merge(grid, biocounts,by.x="cell",by.y="cell")

#Make adjustments so the output is more visually interesting
grid$count    <- log(grid$count)
cutoff        <- quantile(grid$count,0.9)
grid          <- grid %>% mutate(count=ifelse(count>cutoff,cutoff,count))

#Get polygons for each country of the world
countries <- map_data("world")

#Plot everything on a flat map
(p<- ggplot() + 
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid,      aes(x=long, y=lat, group=group, fill=count), alpha=0.4)    +
  geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  geom_point  (aes(x=cellcenters$lon_deg, y=cellcenters$lat_deg)) +
  scale_fill_gradient(low="blue", high="red"))

#Replot on a spherical projection
p+coord_map("ortho", orientation = c(-38.49831, -179.9223, 0))+
  xlab('')+ylab('')+
  theme(axis.ticks.x=element_blank())+
  theme(axis.ticks.y=element_blank())+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_blank())+
  ggtitle('Your data could look like this')

ggsave(p, filename = "outputs/map_cells.png",
       height = 5, width = 8)

# global grid with 100km² ----
#Construct a global grid with cells approximately 1000 miles across
dggs <- dgconstruct(res = 12, metric=FALSE, resround='down')


#Get the corresponding grid cells for each earthquake epicenter (lat-long pair)
bio_short$cell <- dgGEO_to_SEQNUM(dggs,bio_short$LONGITUDE,bio_short$LATITUDE)$seqnum

#Get the number of earthquakes in each equally-sized cell
biocounts   <- bio_short %>% group_by(cell) %>% summarise(count=n())

# joining dataset with all important variables ----
data1 <- bio_full_acc %>% 
  mutate(hpdscale = bio_full_hpd$scalehpd,
         gridcell = bio_short$cell)

# center duration ----

bio_509 <- bio %>% 
  filter(STUDY_ID == "10")

bio_m_509 <- biotime_meta %>% 
  filter(STUDY_ID == "10")

# extracting values online
SP1 <- SP %>% 
  dplyr::select(-STUDY_ID_PLOT) %>% 
  rename( latitude = "LATITUDE",
          longitude = "LONGITUDE")
#write.csv(SP1, "outputs/SP1.csv")
# results in same values

# checking distributions
hist(bio_hpd_short$e_hpd)
max(bio_hpd_short$e_hpd)
min(bio_hpd_short$e_hpd)

hist(log(bio_hpd_scale$scalehpd))
hist(log(bio_aa_scale$scaleacc))
hist(log(bio_aa$e))

# figuring out what to do with biomass ----

unique(bio$ABUNDANCE_TYPE) 
# Count, Density, MeanCount, Presence/Absence
unique(biotime_meta$BIOMASS_TYPE)
# cover, size, volume weight

p_a <- biotime_meta %>% 
  filter(STUDY_ID == "509")

abundance_type_plot <- bio %>% 
  group_by(ABUNDANCE_TYPE) %>% 
  summarise(plots =length(unique(STUDY_ID_PLOT)))

abundance_ab <- bio %>% 
  group_by(AB_BIO) %>% 
  summarise(plots =length(unique(STUDY_ID_PLOT)))

den <- bio %>% 
  filter(AB_BIO  %in% c( "A")) #%>% 
  #filter(sum.allrawdata.ABUNDANCE < 1)

plo <- bio %>% 
  filter(STUDY_ID == "502") %>% 
  distinct(STUDY_ID_PLOT)

den2 <- bio_meta %>% 
  filter(!ABUNDANCE_TYPE == "Count")

# loop with AB = B ----
bio$sum.allrawdata.ABUNDANCE[bio$sum.allrawdata.ABUNDANCE < 1] <- 1

beta_Jacca

# acc testing distribution of fake lat/long ----
library(generator)
fake_lat <- r_latitudes(5800)
fake_long <- r_longitudes(5800)

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

#f_bio_aa_scale <- f_bio_aa_short %>%
#  mutate(f_scaleacc=(f_e-min(f_e))/(max(f_e)-min(f_e)))

hist(f_bio_aa_short$f_e) # more normal distributed?
hist(log(f_bio_aa_scale$f_scaleacc))

# hpd testing distribution of fake lat/long ----
library(generator)
#fake_lat <- r_latitudes(1023)
#fake_long <- r_longitudes(1023)

#fake_ll <- SP %>% 
 # mutate(fake_lat = c(fake_lat),
  #       fake_long = c(fake_long)) %>% 
  #dplyr::select(- LATITUDE, -LONGITUDE)

#f_points <- cbind(fake_ll$fake_long, fake_ll$fake_lat)


#f_sppoints <- SpatialPoints(f_points, proj4string=CRS('+proj=longlat +datum=WGS84'))
f_tp_hpd <- spTransform(f_sppoints, crs(hpd))

f_e_hpd <- extract(hpd, f_tp_hpd)

f_bio_hpd <- cbind(fake_ll, f_e_hpd)
f_bio_hpd_short <- na.omit(f_bio_hpd)

#f_bio_aa_scale <- f_bio_aa_short %>%
#  mutate(f_scaleacc=(f_e-min(f_e))/(max(f_e)-min(f_e)))

hist(f_bio_hpd_short$f_e_hpd) # more normal distributed?
hist(log(f_bio_aa_scale$f_scaleacc))

taxa_st <- data1 %>% 
  group_by(TAXA) %>% 
  summarise(length(unique(STUDY_ID_PLOT)))

check <- data1 %>% 
  filter(scaleacc_25 > 0.5)
rd$Jtu <- as.numeric(beta_Jaccard$Jtu)
hist(beta_Jaccard$Jtu)


