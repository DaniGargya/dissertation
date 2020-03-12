# Analysis
# Dani Gargya, daniela@gargya.de
# March 2020


# calculating jaccard turnover ----
# Load packages ----
library(tidyverse) # (Contains loads of useful functions)
library(ggplot2) # (Package for making nice graphs)

# load data ----
avatar_data <- read.csv("data/avatar_fauna.csv")


# bio <- read.csv("data/bio.csv")
bio_turnover <- bio %>% 
  select(SAMPLE_DESC, STUDY_ID_PLOT, STUDY_ID, PLOT, YEAR, GENUS_SPECIES, ID_SPECIES, sum.allrawdata.ABUNDANCE) %>% 
  group_by(STUDY_ID_PLOT) %>% 
  filter(YEAR %in% c(max(YEAR), min(YEAR))) %>% 
  group_by(STUDY_ID_PLOT, YEAR, GENUS_SPECIES) %>% 
  summarise(Abundance = sum(sum.allrawdata.ABUNDANCE))
  
bio_t_matrix <- bio_turnover %>% 
  group_by(STUDY_ID_PLOT) %>% 
  spread(GENUS_SPECIES, Abundance)







# from biotime website ----
getBeta<-function(x, origin=1, consec=F){
  size<-dim(x)[1]-1
  ifelse(consec==T, id<-cbind(1:size, 1:size+1), id<-cbind(1:dim(x)[1], origin))
  output<-data.frame(Jaccard=as.matrix(1-vegdist(x, method="jaccard"))[id])
  ifelse(consec==T, rownames(output)<-rownames(x)[-1], rownames(output)<-rownames(x))
  output
}

#### get the beta diversity metrics
betaX<-lapply(dsList, function(x)getBeta(x[,-1]))
beta<-do.call(rbind, lapply(betaX, function(x)x[-1,]))
beta$Year<-substring(rownames(beta), nchar(rownames(beta))-3)


