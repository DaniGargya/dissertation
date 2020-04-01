# Question script
# Dani Gargya
# March 2020

# Load libraries ----
library(tidyverse)
library(vegan)
library(betapart)

# workflow ----
  # create species matrix per STUDY_ID_PLOT
  # create binary matrix per STUDY_ID_PLOT
  # calculate jaccard using betapart
  # extract jtu values (turnover component)
  # create dataframe with STUDY_ID_PLOT and matching jtu value

# load data ----
# this is a reduced dataset with my inclusion criteria
# additionally I have summed the abundances for the first and last year each in every plot
# to keep the dataset small for now, I have only included three plots (10_1, 10_2, 10_9)
bio_turnover <- read.csv("data/bio_turnover.csv")

# Remove unnecessarily columns
bio_turnover <- bio_turnover %>% dplyr::select(-X)

# with loop ----

beta_JTU <- data.frame(matrix(ncol = 2, nrow = length(unique(bio_turnover$STUDY_ID_PLOT)))) 
names(beta_JTU) <- c("STUDY_ID_PLOT", "beta_JTU") 
i = 1

# for loop with betapart ----
for (i in 1:length(unique(bio_turnover$STUDY_ID_PLOT))) {
  StudyIDPlot <- as.character(unique(bio_turnover$STUDY_ID_PLOT)[i])
  sub_bio_turnover <- filter(bio_turnover, 
                             STUDY_ID_PLOT == StudyIDPlot) 
  sub_bio_turnover2 <- pivot_wider(sub_bio_turnover, names_from = GENUS_SPECIES, 
                                   values_from = Abundance, 
                                   values_fill = list(Abundance = 0))
  sub_bio_turnover3 <- dplyr::select(sub_bio_turnover2, -STUDY_ID_PLOT, -YEAR) 
  beta_matrix <- with(sub_bio_turnover3, ifelse(sub_bio_turnover3 > 0,1,0))
  beta.JTU <- beta.sample(beta_matrix, index.family="jaccard")$sampled.values$beta.JTU
  beta_JTU[i,] <- c(StudyIDPlot, beta.JTU)
  i = i+1
}
