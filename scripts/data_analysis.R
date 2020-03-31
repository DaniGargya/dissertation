# Data analysis for dissertation
# Dani Gargya
# March 20

# workflow ----
# creating a model
# checking it
# sensitivtiy analysis (extra script?)

# loading libraries ----
library(tidyverse) # (Contains loads of useful functions)
library(ggplot2) # (Package for making nice graphs)
library(MCMCglmm)
library(brms)

# importing data ----
#data1 <- read.csv("")

# checking data ----
head(data1)
tail(data1)
sumamry(data1)

# visualising distributions of data----
# -> zero one inflated beta ditribution

# running MCMCglmm mixed effects model ----

# default priors
# pr = TRUE to do predictions later?
# zero one inflated beta distribution
# which effects have random slope/intercept?

mod1 <- MCMCglmm(temporal_turnover ~ accessibility + duration + (area?) + accessibility:HPD + (accessibility|taxa),
                 random = cell, STUDY_ID,
                 data = data1,
                 pr = TRUE, nitt = 100 000, family = 

# with brms? ----                   
# Turnover before/after

# Set priors for a zero-one inflated beta model
# This model was chosen because turnover is bound between 0 and 1
# There are zeros and ones in the data but also continuous values in between
prior2b <- c(set_prior(prior = 'normal(0,6)', class='b', coef='min_years.scaled'), 	# global slope
             set_prior(prior = 'normal(0,6)', class='Intercept', coef=''), 		# global intercept
             set_prior("normal(0,.5)", class = "Intercept", dpar = "zoi"),
             set_prior("normal(0,.5)", class = "Intercept", dpar = "coi"), 	
             set_prior(prior = 'cauchy(0,2)', class='sd'))	# group-level intercepts

# zoi refers to the probability of being a zero or a one

# coi refers to the conditional probability of being a one 
# (given an observation is a zero or a one)

# phi is the precision parameter of zero-one inflated beta distribution
ba_tu <- brm(bf(Jtu_base ~ period + min_years.scaled + 
                  (1|Biome), coi ~ 1, zoi ~ 1),
             family = zero_one_inflated_beta(), 
             data = beta_periods,
             prior = prior2b, iter = 6000,
             warmup = 2000,
             inits = '0',
             control = list(adapt_delta = 0.85),
             cores = 2, chains = 2)

# Check model and save output ----
summary(ba_tu)
plot(ba_tu)
save(ba_tu, file = "data/output/ba_tu2018.RData")