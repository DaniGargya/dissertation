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
             iter = 2000,
             warmup = 1000,
             inits = '0',
             control = list(adapt_delta = 0.85),
             cores = 2, chains = 2)

# Check model and save output ----
summary(mo_tu)
plot(mo_tu)
save(mo_tu, file = "outputs/mo_tu.RData")

# predicting ----
predictions <- ggpredict(mo_tu, terms = c("scaleacc", "scalehpd", "duration_plot_center"))

#predictions$Monitoring <- factor(predictions$group, levels = c("11.44", "20.79", "30.13"),
                                 #labels = c("Short-term", "Moderate", "Long-term"))

plot(predictions)

# saving model output as tabe
table_m <- tidy_stan(mo_tu, effects = "fixed", digits = 2) %>% print()

# Create table for report
#stargazer(table_m, type = "html", out = "outputs/table_m.html", title = "results")

# visualising ----
theme_clean <- function(){
  theme_bw() +
    theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 14, face = "plain"),             
          axis.title.y = element_text(size = 14, face = "plain"),             
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size = 15, vjust = 1, hjust = 0.5),
          legend.text = element_text(size = 12, face = "italic"),          
          legend.title = element_text(size = 12, face = "bold"),                              
          legend.position = c(0.2, 0.8))
}

ggplot() +
  geom_line(data = predictions, aes(x = x, y = predicted),
            size = 2) +
  geom_ribbon(data = predictions, aes(ymin = conf.low, ymax = conf.high, 
                                      x = x), alpha = 0.1) +
  geom_point(data = data1, aes(x = scaleacc, y = Jtu),
             alpha = 0.1, size = 2) +
  #annotate("text", x = -0.65, y = 5, label = "Slope = -0.06, Std. error = 0.01") +  
  #scale_x_continuous(limits = c (0.8, 1)) +
  theme_clean() +
  #scale_fill_manual(values = c("darksalmon", "firebrick3", "firebrick4")) +
  #scale_colour_manual(values = c("darksalmon", "firebrick3", "firebrick4")) +
  labs(x = "\nAccessibility", y = "Jaccard turnover\n")
