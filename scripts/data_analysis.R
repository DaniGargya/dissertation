# Data analysis for dissertation
# Dani Gargya
# March 20

### It's model 6 that I have een trying to run

# loading libraries ----
library(tidyverse) # (Contains loads of useful functions)
library(brms)
#library(rstan)
#library(rstantools)
#library(ggeffects)
#ibrary(stargazer)  # for tables of model outputs
#library(broom)

#library(tidybayes)
#library(bayesplot)
#library(modelr)
#library(sjstats)

# importing data ----
data1 <- read.csv("data/data1.csv") %>%  dplyr::select(-X)

# checking data ----
str(data1)
head(data1)
tail(data1)
sumamry(data1)

# running models ----

# models to work with ----
# main model
mo_tu_simp6 <- brm(bf(Jtu ~ scaleacc_25 + scalehpd_25 + duration_plot + TAXA +
                       (1|cell) + (1|STUDY_ID)),
                   family = zero_one_inflated_beta(), 
                   data = data1,
                   iter = 4000,
                   warmup = 1000,
                   inits = '0',
                   control = list(adapt_delta = 0.85),
                   cores = 4, chains = 4)

save(mo_tu_simp6, file = "outputs/mo_tu_simp6.RData")

# sensitivity analysis plants
data1_plants <- data1 %>% 
  filter(TAXA == "Terrestrial plants")

mo_tu_simp_plants <- brm(bf(Jtu ~ scaleacc_25 + scalehpd_25 + duration_plot +
                        (1|STUDY_ID)), # or with (1|cell)?
                   family = zero_one_inflated_beta(), 
                   data = data1_plants,
                   iter = 4000,
                   warmup = 1000,
                   inits = '0',
                   control = list(adapt_delta = 0.85),
                   cores = 4, chains = 4)

save(mo_tu_simp_plants, file = "outputs/mo_tu_simp_plants.RData")

# sensitivity analysis scale (1km) # but lots of missing data!
mo_tu_simp8 <- brm(bf(Jtu ~ scaleacc_1 + scalehpd_1 + duration_plot + TAXA +
                        (1|cell) + (1|STUDY_ID)), 
                   family = zero_one_inflated_beta(), 
                   data = data1,
                   iter = 4000,
                   warmup = 1000,
                   inits = '0',
                   control = list(adapt_delta = 0.85),
                   cores = 4, chains = 4)

save(mo_tu_simp8, file = "outputs/mo_tu_simp8.RData")


# sensitivity analysis scale (50km)
# done
mo_tu_simp7 <- brm(bf(Jtu ~ scaleacc_50 + scalehpd_50 + duration_plot + TAXA +
                        (1|STUDY_ID)), # or with (1|cell)?
                   family = zero_one_inflated_beta(), 
                   data = data1,
                   iter = 4000,
                   warmup = 1000,
                   inits = '0',
                   control = list(adapt_delta = 0.85),
                   cores = 4, chains = 4)

save(mo_tu_simp7, file = "outputs/mo_tu_simp7.RData")

# sensitivity analysis scale 100km
mo_tu_scale100 <- brm(bf(Jtu ~ scaleacc_100 + scalehpd_100 + duration_plot + TAXA +
                           (1|cell) + (1|STUDY_ID)),
                      family = zero_one_inflated_beta(), 
                      data = data1,
                      iter = 4000,
                      warmup = 1000,
                      inits = '0',
                      control = list(adapt_delta = 0.85),
                      cores = 2, chains = 4)

save(mo_tu_scale100, file = "outputs/mo_tu_scale100.RData")


# sensitivity analysis smaller time frame 1970 - 2010
data1_tfsmall <- data1 %>% 
  filter(!START_YEAR < 1970) %>% 
  filter(!END_YEAR > 2015)

mo_tu_simp_tf <- brm(bf(Jtu ~ scaleacc_25 + scalehpd_25 + duration_plot +
                          (1|cell) + (1|STUDY_ID)), 
                     family = zero_one_inflated_beta(), 
                     data = data1_tfsmall,
                     iter = 4000,
                     warmup = 1000,
                     inits = '0',
                     control = list(adapt_delta = 0.85),
                     cores = 4, chains = 4)

# richness change
mo_tu_simp_ri <- brm(bf(richness_change ~ scaleacc_25 + scalehpd_25 + duration_plot + TAXA +
                        (1|cell) + (1|STUDY_ID)),
                   family = zero_one_inflated_beta(), 
                   data = data1,
                   iter = 4000,
                   warmup = 1000,
                   inits = '0',
                   control = list(adapt_delta = 0.85),
                   cores = 4, chains = 4)

save(mo_tu_simp_ri, file = "outputs/mo_tu_simp_ri.RData")


# predicting ----
(plot <- data1 %>%
    data_grid(scaleacc_25 = seq_range(scaleacc_25, n = 3), scalehpd_25 = seq_range(scalehpd_25, n = 3), duration_plot = seq_range(duration_plot, n = 3), TAXA = rep(c("Birds", "Mammals", "Terrestrial invertebrates", "Terrestrial plants"), 3)) %>%
    add_predicted_draws(mo_tu_simp6, re_formula = NULL, allow_new_levels = TRUE) %>%
    ggplot(aes(x = scaleacc_25)) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .8, .5), colour = "#578988", alpha = 0.5) +
    geom_hline(linetype = "dashed", yintercept = 0, colour = "grey10") +
    geom_point(aes(y = Jtu), data = data1, colour = "#578988",
               alpha = 0.8, size = 2) +
    scale_fill_manual(values = c("grey90", "grey80", "grey60")) +
    labs(x = "\nAccessibility (proportion)", 
         y = "Turnover\n", title = "Dani's plot\n") +
    #scale_x_continuous(breaks = c(-1.727407, -0.7744444, 0.6537037, 2.081852, 3.510000),
    #                   labels = paste0(c("0", "0.08", "0.16", "0.24", "0.32"))) +
    #scale_y_continuous(breaks = c(0, 0.5, 1),
    #                   labels = c("0", "0.5", "1")) +
    theme_classic() +
    guides(fill = F))


# TAXA
(plot2 <- data1 %>%
    data_grid(scaleacc_25 = seq_range(scaleacc_25, n = 3), scalehpd_25 = seq_range(scalehpd_25, n = 3), duration_plot = seq_range(duration_plot, n = 3), TAXA = rep(c("Birds", "Mammals", "Terrestrial invertebrates", "Terrestrial plants"), 3)) %>%
    add_predicted_draws(mo_tu_simp6, re_formula = NULL, allow_new_levels = TRUE) %>%
    ggplot(aes(x = scaleacc_25)) +
    stat_lineribbon(aes(y = .prediction, color = TAXA), .width = c(.95), alpha = 0.4) +
    geom_hline(linetype = "dashed", yintercept = 0, colour = "grey10") +
    geom_point(aes(y = Jtu), data = data1, colour = "#578988",
               alpha = 0.8, size = 2) +
    scale_fill_manual(values = c("grey90", "grey80", "grey60")) +
    labs(x = "\nAccessibility (proportion)", 
         y = "Turnover\n", title = "Dani's plot\n") +
    #scale_x_continuous(breaks = c(-1.727407, -0.7744444, 0.6537037, 2.081852, 3.510000),
    #                   labels = paste0(c("0", "0.08", "0.16", "0.24", "0.32"))) +
    #scale_y_continuous(breaks = c(0, 0.5, 1),
    #                   labels = c("0", "0.5", "1")) +
    theme_classic() +
    guides(fill = F))


# Code from Gergana:

# the variables need to have the exact same names
# slopes_luh_Jtu are the raw points - the exact dataframe used to run the model
# Jtu_luh_cont is the model output

# (turnover_luh2 <- slopes_luh_Jtu %>%
#    data_grid(forest.diff_scaled = seq_range(forest.diff_scaled, n = 101)) %>%
#    add_predicted_draws(Jtu_luh_cont, re_formula = NULL, allow_new_levels = TRUE) %>%
#    ggplot(aes(x = forest.diff_scaled)) +
#    stat_lineribbon(aes(y = .prediction), .width = c(.95, .8, .5), colour = "#578988", alpha = 0.5) +
#    geom_hline(linetype = "dashed", yintercept = 0, colour = "grey10") +
#    geom_point(aes(y = final_tu), data = slopes_luh_Jtu, colour = "#578988",
#               alpha = 0.8, size = 2) +
#    scale_fill_manual(values = c("grey90", "grey80", "grey60")) +
#    labs(x = "\nForest loss (proportion)", 
#         y = "Turnover\n", title = "LUH (full time series duration)\n") +
#    scale_x_continuous(breaks = c(-1.727407, -0.7744444, 0.6537037, 2.081852, 3.510000),
#                       labels = paste0(c("0", "0.08", "0.16", "0.24", "0.32"))) +
#    scale_y_continuous(breaks = c(0, 0.5, 1),
#                       labels = c("0", "0.5", "1")) +
#    guides(fill = F))

# Gergana's code implemented
load("~/Desktop/mo_tu_simp6.RData")
load("outputs/mo_tu_simp1.RData")


# histogram studyID deviation ----
re_studyid <- as.data.frame(ranef(mo_tu_simp6)$STUDY_ID)
hist(re_studyid$Estimate.Intercept)

re_studyid$STUDY_ID <- rownames(re_studyid)
re_studyid$STUDY_ID <- as.factor(re_studyid$STUDY_ID)


re_studyid_ta <- re_studyid %>% 
  left_join(data1, by = "STUDY_ID")

(re_study_hist <- ggplot(re_studyid_ta, aes(x=Estimate.Intercept, color = TAXA, fill = TAXA))+
    geom_histogram() +
    theme_classic())
ranef(mo_tu_simp6)


#### other analyis ----
# check distribution 0,1, 0-1 raw data ----
dis1 <- data1 %>% 
  filter(Jtu == 1) %>% 
  summarise(time_series = length(unique(STUDY_ID_PLOT)))
# 394

dis0 <- data1 %>% 
  filter(Jtu == 0) %>% 
  summarise(time_series = length(unique(STUDY_ID_PLOT)))
# 2123

# between 0 and 1 -> 3270

# quantile analysis ----
# https://strengejacke.github.io/ggeffects/articles/effectsatvalues.html
quantile(data1$scalehpd_25, probs = c(0.3, 0.5, 0.7))
pred_q30 <- ggpredict(mo_tu_simp1,terms = c("scaleacc_25", "scalehpd_25[0.9943990, 0.9966642, 0.9999365]"))

quantile(data1$scalehpd_25, probs = c(0.2, 0.5, 0.8))
pred_q20 <- ggpredict(mo_tu_simp1,terms = c("scaleacc_25", "scalehpd_25[0.9935632, 0.9966642, 0.9999492]"))

quantile(data1$scalehpd_25, probs = c(0.1, 0.5, 0.9))
pred_q10 <- ggpredict(mo_tu_simp1,terms = c("scaleacc_25", "scalehpd_25[0.9868194, 0.9966642, 0.9999771]"))

quantile(data1$scalehpd_25, probs = c(0.02, 0.5, 0.98))
pred_q2 <- ggpredict(mo_tu_simp1,terms = c("scaleacc_25", "scalehpd_25[0.8859076, 0.9966642, 0.9999998]"))

# q30
(plotq30 <- ggplot() +
  geom_line(data = pred_q30, aes(x = x, y = predicted, color = group),
            size = 1) +
  geom_ribbon(data = pred_q30, aes(ymin = conf.low, ymax = conf.high, 
                                        x = x), alpha = 0.1)) #+
  #facet_wrap("group"))

# q20
(plotq20 <- ggplot() +
    geom_line(data = pred_q20, aes(x = x, y = predicted, color = group),
              size = 1) +
    geom_ribbon(data = pred_q20, aes(ymin = conf.low, ymax = conf.high, 
                                     x = x), alpha = 0.1)) #+
  #facet_wrap("group"))


# q10
(plotq10 <- ggplot() +
    geom_line(data = pred_q10, aes(x = x, y = predicted, color = group),
              size = 1) +
    geom_ribbon(data = pred_q10, aes(ymin = conf.low, ymax = conf.high, 
                                     x = x), alpha = 0.1)) #+
#facet_wrap("group"))

# q2
(plotq2 <- ggplot() +
    geom_line(data = pred_q2, aes(x = x, y = predicted, color = group),
              size = 1) +
    geom_ribbon(data = pred_q2, aes(ymin = conf.low, ymax = conf.high, 
                                     x = x), alpha = 0.1))# +
    #facet_wrap("group"))

grid.arrange(plotq30, plotq20, plotq10, plotq2, ncol = 1)


# models that didn't converge ----
# one fixed effect
mo_tu_simp1 <- brm(bf(Jtu ~ scaleacc_25 + scalehpd_25 + duration_plot + TAXA +
                        (1|STUDY_ID)),
                   family = zero_one_inflated_beta(), 
                   data = data1,
                   iter = 4000,
                   warmup = 1000,
                   inits = '0',
                   control = list(adapt_delta = 0.85),
                   cores = 4, chains = 4)

save(mo_tu_simp1, file = "outputs/mo_tu_simp1.RData")

# interaction
mo_tu_simp2 <- brm(bf(Jtu ~ scaleacc_25*scalehpd_25 + duration_plot + TAXA +
                        (1|STUDY_ID)),
                   family = zero_one_inflated_beta(), 
                   data = data1,
                   iter = 4000,
                   warmup = 1000,
                   inits = '0',
                   control = list(adapt_delta = 0.85),
                   cores = 4, chains = 4)

save(mo_tu_simp2, file = "outputs/mo_tu_simp2.RData")

# only interaction
mo_tu_simp3 <- brm(bf(Jtu ~ scaleacc_25*scalehpd_25 +
                        (1|STUDY_ID)),
                   family = zero_one_inflated_beta(),
                   data = data1,
                   iter = 4000,
                   warmup = 1000,
                   inits = '0',
                   control = list(adapt_delta = 0.85),
                   cores = 4, chains = 4)

save(mo_tu_simp3, file = "outputs/mo_tu_simp3.RData")

# interaction, other fixed effect
mo_tu_simp4 <- brm(bf(Jtu ~ scaleacc_25*scalehpd_25 +
                        (1|cell)),
                   family = zero_one_inflated_beta(),
                   data = data1,
                   iter = 4000,
                   warmup = 1000,
                   inits = '0',
                   control = list(adapt_delta = 0.85),
                   cores = 4, chains = 4)

save(mo_tu_simp4, file = "outputs/mo_tu_simp4.RData")

# random slopes taxa
mo_tu_simp9 <- brm(bf(Jtu ~ scaleacc_25 + scalehpd_25 + duration_plot + 
                        (scaleacc_25|TAXA) + (1|cell) + (1|STUDY_ID)),
                   family = zero_one_inflated_beta(), 
                   data = data1,
                   iter = 4000,
                   warmup = 1000,
                   inits = '0',
                   control = list(adapt_delta = 0.85),
                   cores = 4, chains = 4)

save(mo_tu_simp9, file = "outputs/mo_tu_simp9.RData")

# interaction
mo_tu <- brm(bf(Jtu ~ scaleacc*scalehpd + duration_plot_center +
                  (scaleacc|TAXA) + (1|cell) + (1|STUDY_ID), coi ~ 1, zoi ~ 1),
             family = zero_one_inflated_beta(), 
             data = data1,
             iter = 2000,
             warmup = 1000,
             inits = '0',
             control = list(adapt_delta = 0.85),
             cores = 2, chains = 2)

save(mo_tu, file = "outputs/mo_tu.RData")

mo_tu2 <- brm(bf(Jtu ~ scaleacc_25*scalehpd_25 + duration_plot_center +
                   (scaleacc_25|TAXA) + (1|cell) + (1|STUDY_ID), coi ~ 1, zoi ~ 1),
              family = zero_one_inflated_beta(), 
              data = data1,
              iter = 2000,
              warmup = 1000,
              inits = '0',
              control = list(adapt_delta = 0.85),
              cores = 2, chains = 2)

save(mo_tu2, file = "outputs/mo_tu2.RData")