# Data transormtaion for dissertation clean
# Dani Gargya
# April 2020

# workflow ----
# loading all relevant datasets
# filter biotime dataset with right criteria
# calculating jaccard
# extracting values accessibility
# extracting values hpd
# creating gobal grid cell variable
# combine all variables
# center duration

# libraries ----
library(tidyverse) # (Contains loads of useful functions)
library(ggplot2) # (Package for making nice graphs)
library(vegan)
library(betapart)
library(rgdal) # to read in and save spatial data
library(raster) # to allow creation, reading, manip of raster data
library(labdsv)
library(dggridR)


# load data ----
# biotime data
biotime_full <- read.csv("data/BioTIMEQuery02_04_2018.csv")
biotime_meta <- read.csv("data/BioTIMEMetadata_02_04_2018.csv")
#bio <- read.csv("data/bio.csv") %>%  dplyr::select(-X)

# accessibility data
aa <- raster("data/2015_accessibility_to_cities_v1.0/2015_accessibility_to_cities_v1.0.tif")
# https://malariaatlas.org/explorer/#/

# human population density data
hpd <- raster("data/gpw_v4_population_density_rev11_2015_30_sec.tif")
# https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11

# filter biotime data ----

### data inclusion criteria
# terrestrial realm
# minimum study duration 5 years
# at least 15 studies per taxa
# no more than 5000 plots per study
# at least 2 survey points in each plot
# filter by has plot?

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
  filter(!STUDY_ID == 298)  %>%  # has 147201 entries; set upper limit
  group_by(STUDY_ID_PLOT) %>% 
  filter(HAS_PLOT == "Y") %>% 
  filter(YEAR %in% c(max(YEAR), min(YEAR))) %>% 
  mutate(number_plots = length(unique(YEAR))) %>% 
  filter(number_plots == 2) %>% # have min and max year per plot only
  dplyr::select(STUDY_ID, PLOT, STUDY_ID_PLOT, START_YEAR, END_YEAR, duration, TAXA, LATITUDE, LONGITUDE, YEAR, sum.allrawdata.ABUNDANCE, GENUS_SPECIES, LATITUDE, LONGITUDE, AREA_SQ_KM, NUMBER_OF_SAMPLES, ABUNDANCE_TYPE, AB_BIO, BIOMASS_TYPE, PROTECTED_AREA) %>% 
  ungroup()

write.csv(bio, "data/bio.csv")

# other useful dataset variations ----
bio_short <- bio %>% 
  distinct(STUDY_ID_PLOT, .keep_all = TRUE) %>% 
  dplyr::select(-YEAR, -sum.allrawdata.ABUNDANCE, -GENUS_SPECIES)

# calculating jacccard ----

# data manipulation ----
# change biomass and density data to 1
bio$sum.allrawdata.ABUNDANCE[bio$sum.allrawdata.ABUNDANCE < 1] <- 1

bio_turnover <- bio %>% 
  dplyr::select(STUDY_ID_PLOT, YEAR, GENUS_SPECIES, sum.allrawdata.ABUNDANCE) %>% 
  group_by(STUDY_ID_PLOT, YEAR, GENUS_SPECIES) %>% 
  summarise(Abundance = sum(sum.allrawdata.ABUNDANCE)) %>% 
  ungroup()

# creating empty dataframe ----
beta_Jaccard <- data.frame(matrix(ncol = 6, nrow = length(unique(bio_turnover$STUDY_ID_PLOT)))) 
names(beta_Jaccard) <- c("STUDY_ID_PLOT", "duration_plot", "richness_change", "Jbeta", "Jtu", "Jne") 
i = 1

# for loop with betapart ----
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

#write.csv(beta_Jaccard, "data/beta_Jaccard_df.csv")
betaj <- read.csv("data/beta_Jaccard_df.csv")
beta_Jaccard <- betaj %>% 
  dplyr::select(-X)

# extracting values accessibility ----
# df of all unique lat/long values
SP <- bio_short %>% 
  dplyr::select(LATITUDE, LONGITUDE, STUDY_ID_PLOT) %>% 
  distinct(LATITUDE, .keep_all = TRUE)

points <- cbind(SP$LONGITUDE, SP$LATITUDE)
sppoints <- SpatialPoints(points, proj4string=CRS('+proj=longlat +datum=WGS84'))
tp <- spTransform(sppoints, crs(aa))

# extract numbers at different scales
e <- extract(aa, tp)
e_2 <- extract(aa, tp, buffer = 2000, fun = mean)
e_5 <- extract(aa, tp, buffer = 5000, fun = mean)
e_25 <- extract(aa, tp, buffer = 25000, fun = mean)
e_50 <- extract(aa, tp, buffer = 50000, fun = mean)
e_75 <- extract(aa, tp, buffer = 75000, fun = mean)
e_100 <- extract(aa, tp, buffer = 100000, fun = mean)

bio_aa <- cbind(SP, e)
bio_aa_2 <- cbind(SP, e, e_2, e_5, e_25, e_50, e_75, e_100)

# save dataframe
#write.csv(bio_aa_2, "data/df_aa_scales.csv")
bio_aa_2 <- read.csv("data/df_aa_scales.csv") %>%  dplyr::select(-X)

# drop NA according to scale I am looking at!!
# add scale
bio_aa_short <- bio_aa_2 %>% 
  drop_na(e_25) %>% 
  mutate(scaleacc_1= 1 - ((e_25 -min(e_25))/(max(e_25)-min(e_25))))

min(bio_aa_short$e_25) # 0.004530478
max(bio_aa_short$e_25) # 8348.096 # only 20%

aa_mm <- setMinMax(aa)
minValue(aa_mm) # 0
maxValue(aa_mm) # 41556
hist(aa)

# over-representing highly accessible places but not full spectrum of not accessible places
# but also most places highly accessible?


hist(bio_aa_short$scaleacc_1)
hist(bio_aa_short$e_25)

bio_full_acc <- bio_aa_short %>% 
  right_join(bio_short, by = "LATITUDE") %>% 
  dplyr::select(-STUDY_ID_PLOT.x, -LONGITUDE.y) %>% 
  rename(STUDY_ID_PLOT = STUDY_ID_PLOT.y,
         LONGITUDE = LONGITUDE.x) %>% 
  dplyr::select(STUDY_ID_PLOT, scaleacc_1, e, e_2, e_5, e_25, e_50, e_75, e_100)


# extracting values hpd ----

# turn lat/long values into right CRS format
tp_hpd <- spTransform(sppoints, crs(hpd))

# extract long/lat from raster at different scales
e_hpd <- extract(hpd, tp_hpd)
e_hpd2 <- extract(hpd, tp_hpd, buffer = 2000, fun = mean)
e_hpd5 <- extract(hpd, tp_hpd, buffer = 5000, fun = mean)
e_hpd25 <- extract(hpd, tp_hpd, buffer = 25000, fun = mean)
e_hpd50 <- extract(hpd, tp_hpd, buffer = 50000, fun = mean)
e_hpd75 <- extract(hpd, tp_hpd, buffer = 75000, fun = mean)
e_hpd100 <- extract(hpd, tp_hpd, buffer = 100000, fun = mean)

# bind extracted values to dataframe
bio_hpd <- cbind(SP, e_hpd)
bio_hpd2 <- cbind(SP, e_hpd, e_hpd2, e_hpd5, e_hpd25, e_hpd50, e_hpd75, e_hpd100)

# save dataframe
#write.csv(bio_hpd2, "data/df_hpd_scales.csv")
bio_hpd2 <- read.csv("data/df_hpd_scales.csv") %>%  dplyr::select(-X)

# drop NA according to scale I am looking at!!
# add scale
bio_hpd_short <- bio_hpd2 %>% 
  drop_na(e_hpd25) %>% 
  mutate(scalehpd_1= (e_hpd25-min(e_hpd25))/(max(e_hpd25)-min(e_hpd25)))

bio_hpd2[!complete.cases(bio_hpd2$e_hpd25),]
# 377_300440

min(bio_hpd_short$e_hpd25) # 0
max(bio_hpd_short$e_hpd25) # 10101.26

hpd_mm <- setMinMax(hpd)
minValue(hpd_mm) # 0
maxValue(hpd_mm) # 483318.2
hist(hpd)

#hist(log(bio_hpd_short$scalehpd_25))
hist(bio_hpd_short$scalehpd_1)

bio_full_hpd <- bio_hpd_short %>% 
  right_join(bio_short, by = "LATITUDE") %>% 
  dplyr::select(-STUDY_ID_PLOT.x, -LONGITUDE.y) %>% 
  rename(STUDY_ID_PLOT = STUDY_ID_PLOT.y,
         LONGITUDE = LONGITUDE.x) %>% 
  dplyr::select(STUDY_ID_PLOT, scalehpd_1, e_hpd, e_hpd2, e_hpd5, e_hpd25, e_hpd50, e_hpd75, e_hpd100)


# creating global grid cell variable ----
#Construct a global grid with cells approximately 100km² (res= 12 equivalates	95.97785km² cell area)
dggs <- dgconstruct(res = 12, metric=FALSE, resround='down')

#Get the corresponding grid cells for each lat-long pair
bio_short$cell <- dgGEO_to_SEQNUM(dggs,bio_short$LONGITUDE,bio_short$LATITUDE)$seqnum

#Get the number of time-series in each equally-sized cell
biocounts   <- bio_short %>% group_by(cell) %>% summarise(count=n())


# joining dataset with all important variables ----
data1 <- beta_Jaccard %>% 
  left_join(bio_full_acc, by = "STUDY_ID_PLOT") %>% 
  left_join(bio_full_hpd, by = "STUDY_ID_PLOT") %>% 
  left_join(bio_short, by = "STUDY_ID_PLOT")
  
  
#filter(!STUDY_ID_PLOT == "377_300440") %>% 
#select(STUDY_ID_PLOT, scaleacc_1, scalehpd_1)
#data1_1$STUDY_ID_PLOT <- as.integer(data1_1$STUDY_ID_PLOT)

# center duration and area----
data1$duration_plot_center <- scale(as.numeric(data1$duration_plot), scale = FALSE)
data1$area_center <- scale(data1$AREA_SQ_KM, scale = FALSE)

# checking and fixing structure of data1 ----
str(data1)
data1$richness_change <- as.numeric(data1$richness_change)
data1$duration_plot <- as.numeric(data1$duration_plot)
data1$Jtu <- as.numeric(data1$Jtu)
data1$STUDY_ID <- as.factor(data1$STUDY_ID)
data1$cell <- as.factor(data1$cell)
data1$STUDY_ID_PLOT <- as.integer(data1$STUDY_ID_PLOT)

#write.csv(data1, "data/data1.csv")