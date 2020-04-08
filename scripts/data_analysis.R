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
library(ggeffects)
library(stargazer)  # for tables of model outputs
library(broom)

library(tidybayes)
library(bayesplot)
library(modelr)
library(sjstats)

# importing data ----
#data1 <- read.csv("")

# checking data ----
head(data1)
tail(data1)
sumamry(data1)

# visualising distributions of data----
# -> zero one inflated beta ditribution

# Model with all ----   
# default priors
# zero one inflated beta distribution
# zoi refers to the probability of being a zero or a one
# coi refers to the conditional probability of being a one 
# (given an observation is a zero or a one)
# phi is the precision parameter of zero-one inflated beta distribution

# one value for each taxa?


mo_tu2 <- brm(bf(Jtu ~ scaleacc_25*scalehpd_25 + duration_plot_center +
                  (scaleacc_25|TAXA) + (1|cell) + (1|STUDY_ID), coi ~ 1, zoi ~ 1),
             family = zero_one_inflated_beta(), 
             data = data1,
             iter = 2000,
             warmup = 1000,
             inits = '0',
             control = list(adapt_delta = 0.85),
             cores = 2, chains = 2)

# Check model and save output ----
summary(mo_tu)
plot(mo_tu)
save(mo_tu, file = "outputs/mo_tu.RData")

pairs(mo_tu)

# predicting ----
predictions <- ggpredict(mo_tu, terms = c("scaleacc", "scalehpd", "duration_plot_center"))

predictions$hpd <- factor(predictions$group, levels = c("-0.08", "0.02", "0.13"),
                                 labels = c("Low", "Moderate", "High"))

plot(predictions)

# saving model output as table
table_m <- tidy_stan(mo_tu, effects = "fixed", digits = 2) %>% print()

# Create table for report
#stargazer(table_m, type = "html", out = "outputs/table_m.html", title = "results")

# calculating quartiles for levels hpd ----
dat <- data1 %>% 
  drop_na(scalehpd)

quantile(dat$scalehpd)
#           0%          25%          50%          75%         100% 
#0.000000e+00 0.000000e+00 8.891233e-05 1.794126e-03 1.000000e+00 

# "0", "8.891233e-05", "1.794126e-03"

quantile(predictions$group)
str(predictions)

# model all with area ----
mo_tu_a <- brm(bf(Jtu ~ scaleacc_25*scalehpd_25 + duration_plot_center + area_center +
                  (scaleacc_25|TAXA) + (1|cell) + (1|STUDY_ID), coi ~ 1, zoi ~ 1),
             family = zero_one_inflated_beta(), 
             data = data1,
             iter = 2000,
             warmup = 1000,
             inits = '0',
             control = list(adapt_delta = 0.85),
             cores = 2, chains = 2)

# model accessibility ----
mo_tu_acc <- brm(bf(Jtu ~ scaleacc + duration_plot_center + area_center +
                    (scaleacc|TAXA) + (1|cell) + (1|STUDY_ID), coi ~ 1, zoi ~ 1),
               family = zero_one_inflated_beta(), 
               data = data1,
               iter = 2000,
               warmup = 1000,
               inits = '0',
               control = list(adapt_delta = 0.85),
               cores = 2, chains = 2)

# model hpd ----
mo_tu_hpd <- brm(bf(Jtu ~ scaleacc + duration_plot_center + area_center +
                      (scaleacc|TAXA) + (1|cell) + (1|STUDY_ID), coi ~ 1, zoi ~ 1),
                 family = zero_one_inflated_beta(), 
                 data = data1,
                 iter = 2000,
                 warmup = 1000,
                 inits = '0',
                 control = list(adapt_delta = 0.85),
                 cores = 2, chains = 2)

# model to get individual taxa response ----
mo_tu_taxa <- brm(bf(Jtu ~ scaleacc_25*scalehpd_25  + TAXA + scaleacc_25:TAXA +
                     (1|cell) + (1|STUDY_ID), coi ~ 1, zoi ~ 1),
                 family = zero_one_inflated_beta(), 
                 data = data1,
                 iter = 2000,
                 warmup = 1000,
                 inits = '0',
                 control = list(adapt_delta = 0.85),
                 cores = 2, chains = 2)
