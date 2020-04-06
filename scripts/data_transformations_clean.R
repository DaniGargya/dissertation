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
# bio <- read.csv("data/bio.csv")

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
  dplyr::select(STUDY_ID, PLOT, STUDY_ID_PLOT, START_YEAR, END_YEAR, duration, TAXA, LATITUDE, LONGITUDE, YEAR, sum.allrawdata.ABUNDANCE, GENUS_SPECIES, LATITUDE, LONGITUDE, AREA_SQ_KM, NUMBER_OF_SAMPLES, ABUNDANCE_TYPE) %>% 
  ungroup()

# other useful dataset variations ----
bio_short <- bio %>% 
  distinct(STUDY_ID_PLOT, .keep_all = TRUE) %>% 
  dplyr::select(-YEAR, -sum.allrawdata.ABUNDANCE, -GENUS_SPECIES)

# calculating jacccard ----
# workflow
### summarising abundance per plot per year per species
### creating empty dataframe
### running loop

# data manipulation ----
bio_turnover <- bio %>% 
  dplyr::select(STUDY_ID, STUDY_ID_PLOT, YEAR, GENUS_SPECIES, sum.allrawdata.ABUNDANCE) %>% 
  filter(STUDY_ID == "509") %>% 
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

#write.csv(beta_Jaccard, "outputs/beta_Jaccard_df.csv")
betaj <- read.csv("outputs/beta_Jaccard_df.csv")
# 502, 509, many 0, many NaN???
# 502 only one species, 509 all zeros in sum.all.rawdataAbundance?
# filter for abundance type?

# extracting values accessibility ----
# df of all unique lat/long values
SP <- bio_short %>% 
  dplyr::select(LATITUDE, LONGITUDE, STUDY_ID_PLOT) %>% 
  distinct(LATITUDE, .keep_all = TRUE)

points <- cbind(SP$LONGITUDE, SP$LATITUDE)
sppoints <- SpatialPoints(points, proj4string=CRS('+proj=longlat +datum=WGS84'))
tp <- spTransform(sppoints, crs(aa))

e <- extract(aa, tp)
# e <- extract(aa, tp, buffer = 2000, fun = mean)

bio_aa <- cbind(SP, e)
bio_aa_short <- na.omit(bio_aa)

bio_aa_scale <- bio_aa_short %>%
  mutate(scaleacc= 1 - ((e-min(e))/(max(e)-min(e))))

#hist(log(bio_aa_scale$scaleacc))

bio_full_acc <- bio_aa %>% 
  left_join(bio_aa_scale, by = "STUDY_ID_PLOT") %>% 
  dplyr::select(-LATITUDE.y, -LONGITUDE.y, -e.y) %>% 
  rename(LATITUDE = LATITUDE.x) %>% 
  right_join(bio_short, by = "LATITUDE") %>% 
  dplyr::select(STUDY_ID_PLOT.y, scaleacc) %>% 
  rename(STUDY_ID_PLOT = STUDY_ID_PLOT.y)


# extracting values hpd ----

# turn lat/long values into right CRS format
tp_hpd <- spTransform(sppoints, crs(hpd))

# extract long/lat from raster
e_hpd <- extract(hpd, tp_hpd)

# bind extracted values to dataframe
bio_hpd <- cbind(SP, e_hpd)

# omit NAs
bio_hpd_short <- na.omit(bio_hpd)

# scale world population extracted
bio_hpd_scale <- bio_hpd_short %>%
  mutate(scalehpd=(e_hpd-min(e_hpd))/(max(e_hpd)-min(e_hpd)))

# check histogram of values
#hist(bio_hpd_scale$scalehpd)
#hist(log(bio_hpd_scale$scalehpd))

# join hpd to all studies
bio_full_hpd <- bio_hpd %>% 
  left_join(bio_hpd_scale, by = "STUDY_ID_PLOT") %>% 
  dplyr::select(-LATITUDE.y, -LONGITUDE.y, -e_hpd.y) %>% 
  rename(LATITUDE = LATITUDE.x) %>% 
  right_join(bio_short, by = "LATITUDE") %>% 
  dplyr::select(STUDY_ID_PLOT.y, scalehpd) %>% 
  rename(STUDY_ID_PLOT = STUDY_ID_PLOT.y)


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

# center duration and area----
data1$duration_plot_center <- scale(as.numeric(data1$duration_plot), scale = FALSE)
data1$area_center <- scale(data1$AREA_SQ_KM, scale = FALSE)

# checking and fixing structure of data1 ----
str(data1)
data1$richness_change <- as.numeric(data1$richness_change)
data1$duration_plot <- as.numeric(data1$duration_plot)
data1$Jtu <- as.numeric(data1$Jtu)
data1$STUDY_ID <- as.factor(data1$STUDY_ID )
data1$cell <- as.factor(data1$cell)

