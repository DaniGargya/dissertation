# Question script
# Dani Gargya
# March 2020

# Load libraries ----
library(tidyverse)
library(vegan)
library(betapart)

# load data ----
# this is a reduced dataset with my inclusion criteria
# additionally I have summed the abundances for the first and last year each in every plot
# to keep the dataset small for now, I have only included three plots (10_1, 10_2, 10_9)
bio_turnover <- read.csv("data/bio_turnover.csv")

# without loop ----
comm <- bio_turnover %>% 
  dplyr::select(-X) %>% 
  spread(GENUS_SPECIES, Abundance, fill =0) %>% 
  dplyr::select(-YEAR, -STUDY_ID_PLOT)

comm_binary <- with(comm, ifelse(comm >0,1,0))

j_components <- beta.pair(comm_binary, index.family = "jaccard")

jtu <- as.matrix (j_components$beta.jtu)

# I get the right numbers, there are just one of many in the jtu matrix
# Right numbers are: 10_1: 0.2, 10_2: 0.181, 10_9: 0.461 (I have checked them individually)
# There is a pattern to the occuring of those right numbers which is 2 to the side and 2 down (see screenshot in comment)
#### How can I extract those right numbers?




# with loop ----
# loop does not work as I don't get the numbers I have obtained when calculating jaccard values per plot individually
# it probably has sth to do with the list and how the loop is going through the list
# since I managed to get right numbers with the method above, I thought I would focus on the code above for now.

# create list for each plot ----
bio_t_list <- split(bio_turnover, unique(bio_turnover$STUDY_ID_PLOT))

# for loop with betapart ----
for (i in seq_along(bio_t_list)) {
  comm_l <- bio_t_list[[i]] %>%  # comm_l is smaller than comm
    dplyr::select(-X) %>% 
    spread(GENUS_SPECIES, Abundance, fill = 0) %>% 
    dplyr::select(-STUDY_ID_PLOT, -YEAR)
  
  comm_binary_l <- with(comm_l, ifelse(comm > 0, 1, 0)) 
  
  j_components_l <- beta.pair(comm_binary_l, index.family = "jaccard")
  
  jtu_l <- as.matrix (j_components_l$beta.jtu)
}
