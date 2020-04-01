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

beta_Jaccard <- data.frame(matrix(ncol = 4, nrow = length(unique(bio_turnover$STUDY_ID_PLOT)))) 
names(beta_Jaccard) <- c("STUDY_ID_PLOT", "Jbeta", "Jtu", "Jne") 
i = 1

# for loop with betapart ----
for (i in 1:length(unique(bio_turnover$STUDY_ID_PLOT))) {
  StudyIDPlot <- as.character(unique(bio_turnover$STUDY_ID_PLOT)[i])
  sub_bio_turnover <- filter(bio_turnover, 
                             STUDY_ID_PLOT == StudyIDPlot)
  sub_bio_turnover_wider <- pivot_wider(sub_bio_turnover, names_from = GENUS_SPECIES, 
                                   values_from = Abundance, 
                                   values_fill = list(Abundance = 0))
  sub_bio_turnover_matrix <- dplyr::select(sub_bio_turnover_wider, -STUDY_ID_PLOT, -YEAR) 
  sub_bio_abundance_matrix <- with(sub_bio_turnover_matrix, ifelse(sub_bio_turnover_matrix > 0,1,0))
  J_components <- beta.pair(sub_bio_abundance_matrix, index.family='jaccard')	# distance
  Jbeta <- J_components$beta.jac
  Jtu <- J_components$beta.jtu
  Jne <- J_components$beta.jne
  beta_Jaccard[i,] <- c(StudyIDPlot, Jbeta, Jtu, Jne)

  i = i+1
}
