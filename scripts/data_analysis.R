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
library(brms)

# importing data ----
#data1 <- read.csv("")

# checking data ----
head(data1)
tail(data1)
sumamry(data1)

# visualising distributions of data----
# -> zero one inflated beta ditribution

# Model with brms ----   
# default priors
# zero one inflated beta distribution
# zoi refers to the probability of being a zero or a one
# coi refers to the conditional probability of being a one 
# (given an observation is a zero or a one)
# phi is the precision parameter of zero-one inflated beta distribution


mo_tu <- brm(bf(Jtu ~ scaleacc*scalehpd + duration_plot_center +
                  (scaleacc|TAXA) + (1|cell) + (1|STUDY_ID), coi ~ 1, zoi ~ 1),
             family = zero_one_inflated_beta(), 
             data = data1,
             iter = 6000,
             warmup = 2000,
             inits = '0',
             #control = list(adapt_delta = 0.85),
             cores = 2, chains = 2)

# Check model and save output ----
summary(ba_tu)
plot(ba_tu)
save(ba_tu, file = "/output/mo_tu.RData")