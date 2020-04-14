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
data1 <- read.csv("data/data1.csv") %>%  dplyr::select(-X)

# checking data ----
head(data1)
tail(data1)
sumamry(data1)

# visualising distributions of data----
# -> zero one inflated beta ditribution


# Model with all mo_tu ----   
# default priors
# zero one inflated beta distribution
# zoi refers to the probability of being a zero or a one
# coi refers to the conditional probability of being a one 
# (given an observation is a zero or a one)
# phi is the precision parameter of zero-one inflated beta distribution

# one value for each taxa?


mo_tu <- brm(bf(Jtu ~ scaleacc*scalehpd + duration_plot_center +
                  (scaleacc|TAXA) + (1|cell) + (1|STUDY_ID), coi ~ 1, zoi ~ 1),
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

# checking funnel plot ----
scaleaccc <- data1 %>% 
  dplyr::select(scaleacc_25)
  
plot(scaleaccc$Slope, I(1/scaleaccc$SE), xlim = c(-100,100), ylim = c(0, 60)) # this makes the funnel plot of slope (rate of change in days/year) and precision (1/SE)

# backtransform data outputs ----
# https://vuorre.netlify.com/post/2019/02/18/analyze-analog-scale-ratings-with-zero-one-inflated-beta-models/

posterior_samples(fit, pars = "b_")[,1:4] %>% 
  mutate_at(c("b_phi_Intercept"), exp) %>% 
  mutate_at(vars(-"b_phi_Intercept"), plogis) %>% 
  posterior_summary() %>% 
  as.data.frame() %>% 
  rownames_to_column("Parameter") %>% 
  kable(digits = 2) 





# model 2, same parameters, better numbers mo_tu2 ----
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
summary(mo_tu2)
plot(mo_tu2)
save(mo_tu2, file = "outputs/mo_tu2.RData")

load("outputs/mo_tu2.RData")

# predicting mo_tu2 ----
predictions_2 <- ggpredict(mo_tu2, terms = c("scaleacc_25","scalehpd_25"))

predictions_2$hpd <- factor(predictions_2$group, levels = c("0.91", "0.99", "1.06"),
                          labels = c("Low", "Moderate", "High"))

plot(predictions_2)

rstantools::posterior_predict(mo_tu2)

# plotting posterior distributions ----
posterior <- as.array(mo_tu2)
dim(posterior)
dimnames(posterior)

color_scheme_set("red")
mcmc_intervals(posterior, pars = c("b_scaleacc_25", "b_scalehpd_25", "b_duration_plot_center", "b_scaleacc_25:scalehpd_25")) + theme_clean() +
  #scale_x_continuous(limits = c(0, 0.2), breaks = c(0, 0.05, 0.10, 0.15, 0.20),
                     #labels = c(0, 0.05, 0.10, 0.15, 0.20)) +
  geom_vline(xintercept = 0)

# model diagnostics ----
# https://bayesat.github.io/lund2018/slides/andrey_anikin_slides.pdf
# observed vs predicted
pp = brms::pp_check(mo_tu2)
pp + theme_bw()

# default plot model predictions
brms::marginal_effects(mo_tu2)

# model fit
fit = fitted(
  mo_tu2
)

plot(conditional_effects(mo_tu2), points = TRUE)

pre_mo_tu2 <- predict(mo_tu2)

# checking model convergence ----
# Check model convergence
mcmc_trace(posterior, pars = c("b_scaleacc_25", "b_scalehpd_25", "b_duration_plot_center", "b_scaleacc_25:scalehpd_25"))
mcmc_trace(posterior, pars = c("sd_TAXA__scaleacc_25"))
plot(mo_tu2)


# model 3, more iterations mo_tu3 ----
mo_tu3 <- brm(bf(Jtu ~ scaleacc_25*scalehpd_25 + duration_plot + AREA_SQ_KM +
                   (scaleacc_25|TAXA) + (1|cell) + (1|STUDY_ID), coi ~ 1, zoi ~ 1),
              family = zero_one_inflated_beta(), 
              data = data1,
              iter = 3000,
              warmup = 1000,
              inits = '0',
              control = list(adapt_delta = 0.85),
              cores = 2, chains = 2)
# Check model and save output ----
summary(mo_tu3)
plot(mo_tu3)
save(mo_tu3, file = "outputs/mo_tu3.RData")
load("outputs/mo_tu3.RData")


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
# did not converge
mo_tu_hpd <- brm(bf(Jtu ~ scalehpd_25 + duration_plot_center + area_center +
                      (scalehpd_25|TAXA) + (1|cell) + (1|STUDY_ID), coi ~ 1, zoi ~ 1),
                 family = zero_one_inflated_beta(), 
                 data = data1,
                 iter = 2000,
                 warmup = 1000,
                 inits = '0',
                 control = list(adapt_delta = 0.85),
                 cores = 2, chains = 2)

# Check model mo_tu_hpd and save output ----
summary(mo_tu_hpd)
plot(mo_tu_hpd)
save(mo_tu_hpd, file = "outputs/mo_tu_hpd.RData")

pairs(mo_tu_hpd)

# model to get individual taxa response ----
mo_tu_taxa <- brm(bf(Jtu ~ scaleacc_25*TAXA +
                     (1|cell) + (1|STUDY_ID), coi ~ 1, zoi ~ 1),
                 family = zero_one_inflated_beta(), 
                 data = data1,
                 iter = 2000,
                 warmup = 1000,
                 inits = '0',
                 control = list(adapt_delta = 0.85),
                 cores = 1, chains = 2)
