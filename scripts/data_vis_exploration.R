# Data visualisation exploration ----
# Dani Gargya
# March 2020

# Load data ----
data1 <- read.csv("data/data1.csv") %>%  dplyr::select(-X)

# load simple model outputs
# no cell
#load("outputs/mo_tu_simp1.RData")
#summary(mo_tu_simp1)

# interaction
# load("outputs/mo_tu_simp2.RData")
#summary(mo_tu_simp2)

# interaction
#load("outputs/mo_tu_simp3.RData")
#summary(mo_tu_simp3)

# interaction
#load("outputs/mo_tu_simp4.RData")
#summary(mo_tu_simp4)

# accessibility|taxa
#load("outputs/mo_tu_simp9.RData")
#summary(mo_tu_simp9)

# MODEL TO USE
load("outputs/mo_tu_simp6.RData")
summary(mo_tu_simp6)
plot(mo_tu_simp6)

# sensitivity 1km
load("outputs/mo_tu_simp8.RData")
summary(mo_tu_simp8)

# sensitivity 50km
load("outputs/mo_tu_simp7.RData")
summary(mo_tu_simp7)

# sensitivity 100km
load("outputs/mo_tu_scale100.RData")
summary(mo_tu_scale100)

# sensitivity plants
load("outputs/mo_tu_simp_plants.RData")
summary(mo_tu_simp_plants)

# richness
load("outputs/mo_tu_simp_ri.RData")
summary(mo_tu_simp_ri)


# Load libraries ----
library(tidyverse) # contains dplyr, ggplot, ...
library(maps) # for mapping the flamingo data using coordinates
library(ggthemes) # for data visualisation
library(treemap) # to create treemaps (taxa boxes)
library(RColorBrewer)
library(treemapify) # for area graph
library(ggplot2)
library(ggpubr)
library(mapdata)
library(gridExtra)
library(ggExtra)
library(tidybayes)
library(ggeffects)

# CREATE THE COLOUR PALETTE
display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n = 8, name = "Dark2")
taxa.palette <- c("#D95F02", "#7570B3", "#E7298A", "#E6AB02") 
names(taxa.palette) <- levels(data1$TAXA)



# clean theme ----
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



# creating niwot theme----
theme_niwot <- function(){
  theme_bw() +
    theme(text = element_text(family = "Helvetica Light"),
          axis.text = element_text(size = 16), 
          axis.title = element_text(size = 18),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 18, vjust = 1, hjust = 0),
          legend.text = element_text(size = 12),          
          legend.title = element_blank(),                              
          legend.position = c(0.95, 0.15), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}

#### RQ1: jaccard ~ accessibility ----
# accessibility model predictions ----
predictions_acc <- ggpredict(mo_tu2, terms = c("scaleacc_25"))
pred_simp1_acc <- ggpredict(mo_tu_simp1, terms = c("scaleacc_25", "duration_plot"))
pred_simp1_acc_re <- ggpredict(mo_tu_simp1, terms = c("scaleacc_25", "duration_plot", "STUDY_ID"), type = "re", allow_new_levels = TRUE)
pred_scale100 <- ggpredict(mo_tu_scale100, terms = c("scaleacc_100"))
pred_scale100_d <- ggpredict(mo_tu_scale100, terms = c("scaleacc_100", "duration_plot"))
pred_scalehpd100 <- ggpredict(mo_tu_scale100, terms = c("scalehpd_100"))
pred_simp6_hpd <- ggpredict(mo_tu_simp6, terms = c("scalehpd_25"))
pred_simp6_hpd_d <- ggpredict(mo_tu_simp6, terms = c("scalehpd_25", "duration_plot"))
pred_simp6_acc <- ggpredict(mo_tu_simp6, terms = c("scaleacc_25", "duration_plot"))
pred_ri <- ggpredict(mo_tu_simp_ri, terms = c("scaleacc_25"))

# visualse pred by study ID ----
library(randomcoloR)
n <- 100
palette <- distinctColorPalette(n)
ggpredict(mo_tu_scale100, terms = c("scaleacc_100", "STUDY_ID"), type = "re") %>% plot(colors = palette) # by study ID
ggpredict(mo_tu_simp6, terms = c("scaleacc_25", "cell"), type = "re") %>% plot(colors = palette) # by study ID


# graph only acc----
(graph_acc <- ggplot() +
  geom_line(data = pred_simp6_acc, aes(x = x, y = predicted),
            size = 2) +
  geom_ribbon(data = pred_simp6_acc, aes(ymin = conf.low, ymax = conf.high, 
                                        x = x), alpha = 0.1) +
  geom_point(data = data1, aes(x = scaleacc_25, y = Jtu, colour = scalehpd_25),
             alpha = 0.1, size = 2) +
  #annotate("text", x = -0.65, y = 5, label = "Slope = -0.06, Std. error = 0.01") +  
  #scale_x_continuous(limits = c (0.8, 1)) +
  theme_clean() +
  labs(x = "\nAccessibility", y = "Turnover\n"))

ggsave(graph_acc, filename = "outputs/graph_simp1_acc.png",
       height = 5, width = 8)

# graph only acc raw----
data_npa <- data1 %>% 
  filter(PROTECTED_AREA == TRUE)

(graph_acc <- ggplot() +
   geom_point(data = data_npa, aes(x = scaleacc_25, y = Jtu, colour = PROTECTED_AREA),
              alpha = 0.2, size = 2) +
   theme_clean() +
   labs(x = "\nAccessibility", y = "Turnover\n"))

(graph_acc <- ggplot() +
    geom_point(data = data1, aes(x = scaleacc_25, y = Jtu, colour = PROTECTED_AREA),
               alpha = 0.5, size = 2) +
    scale_colour_manual(values = c("yellow", "black")) +
    theme_clean() +
    labs(x = "\nAccessibility", y = "Turnover\n"))


# graph only hpd----
(graph_acc <- ggplot() +
   geom_line(data = pred_simp6_hpd, aes(x = x, y = predicted),
             size = 2) +
   geom_ribbon(data = pred_simp6_hpd, aes(ymin = conf.low, ymax = conf.high, 
                                          x = x), alpha = 0.1) +
   geom_point(data = data1, aes(x = scalehpd_25, y = Jtu),
              alpha = 0.5, size = 2) +
   #annotate("text", x = -0.65, y = 5, label = "Slope = -0.06, Std. error = 0.01") +  
   #scale_x_continuous(limits = c (0.8, 1)) +
   theme_clean() +
   labs(x = "\nHPD", y = "Turnover\n"))


# acc and duration ----
ggplot() +
  geom_line(data = pred_simp6_acc, aes(x = x, y = predicted, color = group),
            size = 2) +
  geom_ribbon(data = pred_simp6_acc, aes(ymin = conf.low, ymax = conf.high, 
                                        x = x, fill = group), alpha = 0.1) +
  geom_point(data = data1, aes(x = scaleacc_25, y = Jtu),
             alpha = 0.1, size = 2) +
  #annotate("text", x = -0.65, y = 5, label = "Slope = -0.06, Std. error = 0.01") +  
  #scale_x_continuous(limits = c (0.8, 1)) +
  theme_clean() +
  scale_fill_manual(values = c("darksalmon", "firebrick3", "firebrick4")) +
  scale_colour_manual(values = c("darksalmon", "firebrick3", "firebrick4")) +
  labs(x = "\nAccessibility", y = "Turnover\n")



# (only accessibility, facet_wrap taxa) ----
ggplot() +
    geom_line(data = predictions_2, aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = predictions_2, aes(ymin = conf.low, ymax = conf.high, 
                                          x = x), alpha = 0.1) +
    geom_point(data = data1, aes(x = scaleacc_25, y = Jtu),
               alpha = 0.1, size = 2) +
    facet_wrap("TAXA") +
    #annotate("text", x = -0.65, y = 5, label = "Slope = -0.06, Std. error = 0.01") +  
    #scale_x_continuous(limits = c (0.8, 1)) +
    theme_clean() +
    #scale_fill_manual(values = c("darksalmon", "firebrick3", "firebrick4")) +
    #scale_colour_manual(values = c("darksalmon", "firebrick3", "firebrick4")) +
    labs(x = "\nAccessibility", y = "Jaccard dissimilarity\n")


# graph only hpd----
(graph_acc <- ggplot() +
   geom_line(data = pred_simp6_hpd, aes(x = x, y = predicted),
             size = 2) +
   geom_ribbon(data = pred_simp6_hpd, aes(ymin = conf.low, ymax = conf.high, 
                                          x = x), alpha = 0.1) +
   geom_point(data = data1, aes(x = scalehpd_25, y = Jtu),
              alpha = 0.5, size = 2) +
   #annotate("text", x = -0.65, y = 5, label = "Slope = -0.06, Std. error = 0.01") +  
   #scale_x_continuous(limits = c (0.8, 1)) +
   theme_clean() +
   labs(x = "\nHPD", y = "Turnover\n"))

# graph hpd and duration ----
(graph_acc <- ggplot() +
   geom_line(data = pred_simp6_hpd_d, aes(x = x, y = predicted, color = group),
             size = 2) +
   geom_ribbon(data = pred_simp6_hpd_d, aes(ymin = conf.low, ymax = conf.high, 
                                          x = x, fill = group), alpha = 0.1) +
   geom_point(data = data1, aes(x = scalehpd_25, y = Jtu),
              alpha = 0.5, size = 2) +
   #annotate("text", x = -0.65, y = 5, label = "Slope = -0.06, Std. error = 0.01") +  
   #scale_x_continuous(limits = c (0.8, 1)) +
   theme_clean() +
   labs(x = "\nHPD", y = "Turnover\n"))



# prediction
pred_acc <- ggpredict(mo_tu_simp6, terms = c("scaleacc_25"))

pred_acc_d$Monitoring <- factor(pred_acc_d$group, levels = c("6", "19.1", "32.2"),
                                labels = c("Short-term", "Moderate", "Long-term"))

# model vis
(graph_acc <- ggplot() +
  geom_line(data = pred_acc, aes(x = x, y = predicted),
            size = 2) +
  geom_ribbon(data = pred_acc, aes(ymin = conf.low, ymax = conf.high, 
                                     x = x), alpha = 0.1) +
  geom_point(data = data1, aes(x = scaleacc_25, y = Jtu),
             alpha = 0.1, size = 2) +
  theme_clean() +
  scale_fill_manual(values = c("#1B9E77", "#66A61E", "#7CFC00")) +
  scale_colour_manual(values = c("#1B9E77", "#66A61E", "#7CFC00")) +
  labs(x = "\nAccessibility (proportion)", y = "Turnover\n") +
  theme(legend.position = c(0.2, 0.8)))


#### RQ2: taxa making raincloud plot ----
# basic exploration ----
# boxplot
boxplot(Jtu ~ TAXA, data = data1)

# on scatterplot
(colour_plot <- ggplot(data1, aes(x = scaleacc_100, y = Jtu, colour = TAXA)) +
    geom_point(size = 2) +
    theme_classic())
# -> plant the only one with variation

# facet wrap
(split_plot <- ggplot(aes(x = scaleacc_100, y = Jtu), data = data1) + 
    geom_point() + 
    facet_wrap(~ TAXA))

# Plot the predictions 
(ggplot(pred_scale100) + 
    geom_line(aes(x = x, y = predicted)) +          # slope
    geom_ribbon(aes(x = x, ymin = conf.low , ymax = conf.high), 
                fill = "lightgrey", alpha = 0.5) +  # error band
    geom_point(data = data1,                      # adding the raw data (scaled values)
               aes(x = scaleacc_100, y = Jtu, colour = TAXA)) + 
    labs(x = "acc (scale 100)d)", y = "Turnover") +
    theme_minimal()
)

# call function by Ben Marwick ----
# This code loads the function in the working environment
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

# visualising
(violin_taxa <- ggplot(data1, aes(x = TAXA, y = Jtu)) +
    geom_violin())

# raincloud plot taxa jtu----
(raincloud_taxa <- 
    ggplot(data = data1, 
           aes(x = reorder(TAXA, desc(Jtu)), y = Jtu, fill = TAXA)) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.8) +
    geom_point(aes(y = Jtu, color = TAXA), 
               position = position_jitter(width = 0.15), size = 1, alpha = 0.1) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
    labs(y = "\nJaccard dissimilarity index", x = NULL) +
    guides(fill = FALSE, color = FALSE) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C")) +
    coord_flip() +
    theme_niwot())

ggsave(raincloud_taxa, filename = "outputs/raincloud_taxa_graph.png",
       height = 5, width = 8)

# same, no flip ----
(raincloud_taxa <- 
   ggplot(data = data1, 
          aes(x = reorder(TAXA, desc(Jtu)), y = Jtu, fill = TAXA)) +
   geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.8) +
   geom_point(aes(y = Jtu, color = TAXA), 
              position = position_jitter(width = 0.15), size = 1, alpha = 0.1) +
   geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
   labs(y = "\nJaccard dissimilarity index", x = NULL) +
   guides(fill = FALSE, color = FALSE) +
   scale_y_continuous(limits = c(0, 1)) +
   scale_fill_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C")) +
   scale_colour_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C")) +
   #coord_flip() +
   theme_niwot())
   
# raincloud plot taxa accessibility ----
(raincloud_acc <- 
    ggplot(data = data1, 
           aes(x = reorder(TAXA, desc(scaleacc_25)), y = scaleacc_25, fill = TAXA)) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.8) +
    geom_point(aes(y = scaleacc_25, color = TAXA), 
               position = position_jitter(width = 0.15), size = 1, alpha = 0.1) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
    labs(y = "\nAccessibility score", x = NULL) +
    guides(fill = FALSE, color = FALSE) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C")) +
    scale_colour_manual(values = c("#5A4A6F", "#E47250",  "#EBB261", "#9D5A6C")) +
    coord_flip() +
    theme_niwot())

ggsave(raincloud_acc, filename = "outputs/raincloud_taxa_acc_graph.png",
       height = 5, width = 8)


# raincloud plot taxa hpd ----

# predict acc + taxa ----
# predict taxa ----
me <- ggpredict(mo_tu_scale100, terms = c("scaleacc_100", "TAXA"))
me_6 <- ggpredict(mo_tu_simp6, terms = c("scaleacc_25", "TAXA"))
pred_taxa <-ggpredict(mo_tu_simp6, terms = c("TAXA"))

# plot
(graph_acc <- ggplot() +
   geom_line(data = me_6, aes(x = x, y = predicted, color = group),
             size = 2) +
   geom_ribbon(data = me_6, aes(ymin = conf.low, ymax = conf.high, 
                                           x = x, fill = group), alpha = 0.1) +
   geom_point(data = data1, aes(x = scaleacc_25, y = Jtu, color = TAXA),
              alpha = 0.1, size = 2) +
   #facet_wrap ("TAXA") +
   #annotate("text", x = -0.65, y = 5, label = "Slope = -0.06, Std. error = 0.01") +  
   #scale_x_continuous(limits = c (0.8, 1)) +
   theme_clean() +
   labs(x = "\nAccessibility", y = "Jaccard dissimilarity\n"))

plot(me_6)

# effect size graph ----
# points plot
(ef_taxa_point <- ggplot(taxa, aes(x = TAXA, y = Estimate.scaleacc_25, group = TAXA, color = TAXA)) +
   geom_pointrange(aes(ymin = Estimate.scaleacc_25-Est.Error.scaleacc_25, ymax =  Estimate.scaleacc_25+Est.Error.scaleacc_25)) +
   theme_clean() +
   theme(legend.position = "none"))

# barplot
(ef_taxa <- ggplot(taxa, aes(x = TAXA, y = Estimate.scaleacc_25, fill = TAXA)) +
  geom_bar(stat="identity")+
   geom_errorbar(aes(ymin = Estimate.scaleacc_25-Est.Error.scaleacc_25, ymax =  Estimate.scaleacc_25+Est.Error.scaleacc_25)) +
   theme_clean())

# effect size points plot mo simple1 ----
summary(mo_tu_scale100)
fixef(mo_tu_simp1)
# create dataframe with effectsizes
# the intercept is the predicted mean of the response when all other predictors have the value “0”
tax <- as.data.frame(fixef(mo_tu_scale100, pars = c("Intercept", "TAXAMammals", "TAXATerrestrialinvertebrates", "TAXATerrestrialplants")))
tax$TAXA <- rownames(tax)

me_taxa <- ggpredict(mo_tu_scale100, terms = c("TAXA"))

(ef_taxa_mo_simp <- ggplot(pred_taxa, aes(x = x, y = predicted, color = x)) +
    geom_point(aes(size = 2)) +
    geom_pointrange(aes(ymin = conf.low, ymax =  conf.high), size = 2) +
    theme_clean() +
    theme(legend.position = "none") +
    geom_hline(yintercept=0, linetype="dashed", size=1) +
    coord_flip())



#### RQ 3: interaction hpd ----
# only accessibility, facet_wrap hpd ----
ggplot() +
  geom_line(data = predictions_3, aes(x = x, y = predicted),
            size = 2) +
  geom_ribbon(data = predictions_3, aes(ymin = conf.low, ymax = conf.high, 
                                        x = x), alpha = 0.1) +
  facet_wrap("hpd") +
  geom_point(data = data1, aes(x = scaleacc_25, y = Jtu),
             alpha = 0.1, size = 2) +
  facet_wrap("hpd_q") +
  theme_clean() +
  labs(x = "\nAccessibility", y = "Jaccard dissimilarity\n")

# only accessibility, colour continuous hpd ----
ggplot() +
  geom_line(data = pred_simp6_acc, aes(x = x, y = predicted),
            size = 2) +
  geom_ribbon(data = pred_simp6_acc, aes(ymin = conf.low, ymax = conf.high, 
                                        x = x), alpha = 0.1) +
  geom_point(data = data1, aes(x = scaleacc_25, y = Jtu, colour = scalehpd_25),
             alpha = 0.1, size = 2) +
  #annotate("text", x = -0.65, y = 5, label = "Slope = -0.06, Std. error = 0.01") +  
  #scale_x_continuous(limits = c (0.8, 1)) +
  theme_clean() +
  #scale_fill_manual(values = c("darksalmon", "firebrick3", "firebrick4")) +
  #scale_colour_manual(values = c("darksalmon", "firebrick3", "firebrick4")) +
  labs(x = "\nAccessibility", y = "Jaccard dissimilarity\n")

# make multi panel plot with filtered data sets ----
# filter dataset predicitons ----
pre_low <- pred_scale100_acc_hpd %>% 
  filter(hpd == "Low")

pre_mod <- pred_scale100_acc_hpd %>% 
  filter(hpd == "Moderate")

pre_high <- pred_scale100_acc_hpd %>% 
  filter(hpd == "High")



# filter dataset all ----
data1$hpd_q <- cut(data1$scalehpd_100, 
                   breaks = c(-Inf, 0.2, 0.8, Inf),
                   labels = c("Low", "Moderate", "High"))
dat_low <- data1 %>% 
  filter(hpd_q == "Low")

dat_mod <- data1 %>% 
  filter(hpd_q == "Moderate")

dat_high <- data1 %>% 
  filter(hpd_q == "High")



# make plots ----
# low
(plot_low <- ggplot() +
  geom_line(data = pre_low, aes(x = x, y = predicted),
            size = 2) +
  geom_ribbon(data = pre_low, aes(ymin = conf.low, ymax = conf.high, 
                                        x = x), alpha = 0.1) +
  geom_point(data = dat_low, aes(x = scaleacc_25, y = Jtu, color = TAXA),
             alpha = 0.5, size = 2) +
  theme_clean() +
  theme(legend.position="left") +
  labs(x = "\nAccessibility", y = "Jaccard dissimilarity\n"))

(p_low <- ggMarginal(plot_low, type="density", size = 3, fill = "slateblue"))

# moderate
(plot_mod <- ggplot() +
  geom_line(data = pre_mod, aes(x = x, y = predicted),
            size = 2) +
  geom_ribbon(data = pre_mod, aes(ymin = conf.low, ymax = conf.high, 
                                  x = x), alpha = 0.1) +
  geom_point(data = dat_mod, aes(x = scaleacc_25, y = Jtu, color = TAXA),
             alpha = 0.5, size = 2) +
  theme_clean() +
  theme(legend.position = "none") +
  labs(x = "\nAccessibility", y = "Jaccard dissimilarity\n"))

(p_mod <- ggMarginal(plot_mod, type="density", size = 3, fill = "slateblue"))

# high
(plot_high <- ggplot() +
  geom_line(data = pre_high, aes(x = x, y = predicted),
            size = 2) +
  geom_ribbon(data = pre_high, aes(ymin = conf.low, ymax = conf.high, 
                                  x = x), alpha = 0.1) +
  geom_point(data = dat_high, aes(x = scaleacc_25, y = Jtu, color = TAXA),
             alpha = 0.5, size = 2) +
  theme_clean() +
  theme(legend.position = "none") +
  labs(x = "\nAccessibility", y = "Jaccard dissimilarity\n"))

(p_high <- ggMarginal(plot_high, type="density", size = 3, fill = "slateblue"))

# panel all ----
panel_hpd_acc <- ggarrange(p_low, p_mod, p_high, ncol = 3, align = c("h"))

ggsave(panel_hpd_acc, filename ="outputs/panel_acc_hpd.png",
       height = 4, width = 9)

# marginal disrtibution around ----
# classic plot :
(p <- ggplot(data1, aes(x = scaleacc_25, y= Jtu, colour = TAXA)) +
   geom_point(alpha = 0.5) +
   #scale_size(range = c(0, 1), name="Jaccard turnover") +
   theme_clean() +
   theme(legend.position='right'))
#theme(legend.position="none") 

# marginal density
(p2 <- ggMarginal(p, type="densigram", size = 1, fill = "slateblue"))

# facet wrap and marginal density ----
acc_hpd <- ggplot() +
  geom_line(data = predictions_3, aes(x = x, y = predicted),
            size = 2) +
  geom_ribbon(data = predictions_3, aes(ymin = conf.low, ymax = conf.high, 
                                        x = x), alpha = 0.1) +
  facet_wrap("hpd") +
  geom_point(data = data1, aes(x = scaleacc_25, y = Jtu),
             alpha = 0.1, size = 2) +
  facet_wrap("hpd") +
  theme_clean() +
  labs(x = "\nAccessibility", y = "Jaccard dissimilarity\n") +
  theme(legend.position = "none")

acc_hpd

(p3 <- ggMarginal(acc_hpd, data1, x = scaleacc_25, y = Jtu, type="density", size = 2, fill = "slateblue"))


## acc and hpd, simp1 ----
pred_simp1_acc_hpd <- ggpredict(mo_tu_simp1, terms = c("scaleacc_25", "scalehpd_25[quart]"))
pred_scale100_acc_hpd <- ggpredict(mo_tu_scale100, terms = c("scaleacc_100", "scalehpd_100[0.001818158, 0.009887737, 0.023398489]")) # quantiles at 0.2, 0.5, 0.8
pred_scale100_acc_hpd <- ggpredict(mo_tu_scale100, terms = c("scaleacc_100", "scalehpd_100[0.2, 0.5, 0.8]")) # values at 0.2, 0.5, 0.8
pred_scale100_acc_hpd$hpd <- factor(pred_scale100_acc_hpd$group, levels = c("0.2", "0.5", "0.8"),
                            labels = c("Low", "Moderate", "High"))

ggplot() +
  geom_line(data = pred_scale100_acc_hpd, aes(x = x, y = predicted, color = hpd),
            size = 2) +
  geom_ribbon(data = pred_scale100_acc_hpd, aes(ymin = conf.low, ymax = conf.high, 
                                         x = x, fill = group), alpha = 0.1) +
  #facet_wrap("group") +
  geom_point(data = data1, aes(x = scaleacc_100, y = Jtu),
             alpha = 0.1, size = 2) +
  #annotate("text", x = -0.65, y = 5, label = "Slope = -0.06, Std. error = 0.01") +  
  #scale_x_continuous(limits = c (0.8, 1)) +
  theme_clean() +
  scale_fill_manual(values = c("darksalmon", "firebrick3", "firebrick4")) +
  scale_colour_manual(values = c("darksalmon", "firebrick3", "firebrick4")) +
  labs(x = "\nAccessibility", y = "Turnover\n")


# PCA ----
library(ape)

# clean dataframe
data_small <- data1 %>% 
  dplyr::select(STUDY_ID_PLOT, TAXA, Jtu, scaleacc_25, scalehpd_25, richness_change, AREA_SQ_KM, duration_plot)

data_small <- column_to_rownames(data_small, "STUDY_ID_PLOT")

# colour by taxa ----
res.pca <- PCA(data_small[,-1], scale.unit = TRUE, graph = TRUE)


fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = data_small$TAXA, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "#9400D3", "#228B22"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)

# biplot
fviz_pca_biplot(res.pca, 
                col.ind = data_small$TAXA, palette = "jco", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Taxa") 

# try website ----
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/#pca-data-format

library("FactoMineR")
library(factoextra)
res.pca <- PCA(data_small[,-1], scale.unit = TRUE, graph = TRUE)

print(res.pca)

# eigenvalues/ variances
eig.val <- get_eigenvalue(res.pca)
eig.val

# scree plot
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

# saving results
var <- get_pca_var(res.pca)
var

# Coordinates
head(var$coord)
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)

# correlation cirlce
# Coordinates of variables
head(var$coord, 4)

# variable correlation plot
fviz_pca_var(res.pca, col.var = "black")

# quality of representation
head(var$cos2, 4)

library("corrplot")
corrplot(var$cos2, is.corr=FALSE)

# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(res.pca, choice = "var", axes = 1:2)

# Color by cos2 values: quality on the factor map
fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

# Change the transparency by cos2 values
fviz_pca_var(res.pca, alpha.var = "cos2")

# contributions of variables to PCs
head(var$contrib, 4)
corrplot(var$contrib, is.corr=FALSE)   

# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)

fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 10)

fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

# Change the transparency by contrib values
fviz_pca_var(res.pca, alpha.var = "contrib")

# Color by a custom continuous variable
# Create a random continuous variable of length 10
set.seed(123)
my.cont.var <- rnorm(6)
# Color variables by the continuous variable
fviz_pca_var(res.pca, col.var = my.cont.var,
             gradient.cols = c("blue", "yellow", "red"),
             legend.title = "Cont.Var")

# Color by groups
# Create a grouping variable using kmeans
# Create 3 groups of variables (centers = 3)
set.seed(123)
res.km <- kmeans(var$coord, centers = 3, nstart = 25)
grp <- as.factor(res.km$cluster)
# Color variables by groups
fviz_pca_var(res.pca, col.var = grp, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),
             legend.title = "Cluster")


# dimension description
res.desc <- dimdesc(res.pca, axes = c(1,2), proba = 0.05)
# Description of dimension 1
res.desc$Dim.1

# Description of dimension 2
res.desc$Dim.2

# graph of individuals
ind <- get_pca_ind(res.pca)
ind

# Coordinates of individuals
head(ind$coord)
# Quality of individuals
head(ind$cos2)
# Contributions of individuals
head(ind$contrib)

# plots quality and contribution
fviz_pca_ind(res.pca)

# colour by cos2 values
fviz_pca_ind(res.pca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

# add predicted draws ----
library(modelr)
library(bayesplot)
library(tidybayes)

(model_fit <- data1 %>%
    data_grid(scaleacc_25 = seq_range(scaleacc_25, n = 101)) %>%
    add_predicted_draws(mo_tu_simp1, re_formula = NULL, allow_new_levels = TRUE) %>%
    ggplot(aes(x = scaleacc_25, y = Jtu)) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),
                    alpha = 1/2, colour = "black") +
    geom_point(data = data1, colour = "darkseagreen4", size = 3) +
    scale_fill_brewer(palette = "Greys"))

#mo_tu_simp1_acc <- fixef(mo_tu_simp1, pars = c("scaleacc_25"), summary = FALSE)
summary(mo_tu_simp1)


(plot <- data1 %>%
    data_grid(scaleacc_25 = seq_range(scaleacc_25, n = 3), scalehpd_25 = seq_range(scalehpd_25, n = 3), duration_plot = seq_range(duration_plot, n = 3), TAXA = rep(c("Birds", "Mammals", "Terrestrial invertebrates", "Terrestrial plants"), 3)) %>%
    add_predicted_draws(mo_tu_simp6, re_formula = NULL, allow_new_levels = TRUE) %>%
    ggplot(aes(x = scaleacc_25)) +
    stat_lineribbon(aes(y = .prediction, colour = TAXA, fill = TAXA), .width = c(.95), alpha = 0.5) +
    geom_hline(linetype = "dashed", yintercept = 0, colour = "grey10") +
    geom_point(aes(y = Jtu), data = data1, colour = "#578988",
               alpha = 0.8, size = 2) +
    scale_fill_manual(values = c("grey90", "grey80", "grey60", "darkgrey")) +
    labs(x = "\nAccessibility (proportion)", 
         y = "Turnover\n", title = "Dani's plot\n") +
    #scale_x_continuous(breaks = c(-1.727407, -0.7744444, 0.6537037, 2.081852, 3.510000),
    #                   labels = paste0(c("0", "0.08", "0.16", "0.24", "0.32"))) +
    #scale_y_continuous(breaks = c(0, 0.5, 1),
    #                   labels = c("0", "0.5", "1")) +
    theme_classic() +
    guides(fill = F))

# sensitivity analysis ----
# richness ----
# richness and acc 
(graph_acc <- ggplot() +
   geom_line(data = pred_ri, aes(x = x, y = predicted),
             size = 2) +
   geom_ribbon(data = pred_ri, aes(ymin = conf.low, ymax = conf.high, 
                                          x = x), alpha = 0.1) +
   geom_point(data = data1, colour = "#578988", aes(x = scaleacc_25, y = richness_change),
              alpha = 0.1, size = 2) +
   #annotate("text", x = -0.65, y = 5, label = "Slope = -0.06, Std. error = 0.01") +  
   #scale_x_continuous(limits = c (0.8, 1)) +
   theme_clean() +
   labs(x = "\nAccessibility", y = "Richness change\n"))

# richness and hpd
pred_ri <- ggpredict(mo_tu_simp_ri, terms = c("scalehpd_25"))
(graph_acc <- ggplot() +
    geom_line(data = pred_ri, aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_ri, aes(ymin = conf.low, ymax = conf.high, 
                                    x = x), alpha = 0.1) +
    geom_point(data = data1, colour = "#578988", aes(x = scalehpd_25, y = richness_change),
               alpha = 0.1, size = 2) +
    #annotate("text", x = -0.65, y = 5, label = "Slope = -0.06, Std. error = 0.01") +  
    #scale_x_continuous(limits = c (0.8, 1)) +
    theme_clean() +
    labs(x = "\nHPD", y = "Richness change\n"))

# plants vs all ----
pred_plants <- ggpredict(mo_tu_simp_plants, terms = c("scaleacc_25"))
pred_all <- ggpredict(mo_tu_simp6, terms = c("scaleacc_25"))

data1_plants <- data1 %>% 
  filter(TAXA == "Terrestrial plants")

data1_nopl <- data1 %>% 
  filter(!TAXA == "Terrestrial plants")

(graph_acc <- ggplot() +
    geom_line(data = pred_plants, color = "yellow", aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_plants, color = "yellow", aes(ymin = conf.low, ymax = conf.high, 
                                    x = x), alpha = 0.1, fill = "yellow") +
    geom_point(data = data1_plants, color = "#E7298A", aes(x = scaleacc_25, y = Jtu),
               alpha = 0.1, size = 2) +
    geom_line(data = pred_all, color = "#1B9E77", aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_all, color = "#1B9E77", aes(ymin = conf.low, ymax = conf.high, 
                                        x = x), alpha = 0.1, fill = "#1B9E77") +
    geom_point(data = data1_nopl, color = "#1B9E77", aes(x = scaleacc_25, y = Jtu),
               alpha = 0.1, size = 2) +
    theme_clean() +
    labs(x = "\nAccessibility", y = "Turnover\n"))

# only plants color by study ID ----
(plants_studyID <- ggplot() +
   geom_line(data = pred_plants, color = "yellow", aes(x = x, y = predicted),
             size = 2) +
   geom_ribbon(data = pred_plants, color = "yellow", aes(ymin = conf.low, ymax = conf.high, 
                                                         x = x), alpha = 0.1, fill = "yellow") +
   geom_point(data = data1_plants, aes(x = scaleacc_25, y = Jtu, color = STUDY_ID),
              alpha = 0.5, size = 1) +
   theme_clean() +
   labs(x = "\nAccessibility", y = "Turnover\n") +
   theme(legend.position = "right"))

ggsave(plants_studyID, file = "outputs/plants_studyID.png", width = 7, height = 5)

# scale analysis ----
pred_1 <- ggpredict(mo_tu_simp8, terms = c("scaleacc_1"))
pred_25 <- ggpredict(mo_tu_simp6, terms = c("scaleacc_25"))
pred_50 <- ggpredict(mo_tu_simp7, terms = c("scaleacc_50"))
pred_100 <- ggpredict(mo_tu_scale100, terms = c("scaleacc_100"))

(graph_scale <- ggplot() +
    geom_line(data = pred_1, color = "red", aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_1, color = "red", aes(ymin = conf.low, ymax = conf.high, 
                                                           x = x), alpha = 0.1, fill = "red") +
    geom_line(data = pred_25, color = "green", aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_25, color = "green", aes(ymin = conf.low, ymax = conf.high, 
                                                        x = x), alpha = 0.1, fill = "green") +
    geom_line(data = pred_50, color = "blue", aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_50, color = "blue", aes(ymin = conf.low, ymax = conf.high, 
                                                     x = x), alpha = 0.1, fill = "blue") +
    geom_line(data = pred_100, color = "yellow", aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_100, color = "yellow", aes(ymin = conf.low, ymax = conf.high, 
                                                     x = x), alpha = 0.1, fill = "yellow") +
    theme_clean() +
    labs(x = "\nAccessibility", y = "Turnover\n"))

ggsave(graph_scale, filename = "outputs/graph_scale.png", height = 5, width = 8)

# small time frame analysis ----
pred_full <- ggpredict(mo_tu_simp6, terms = c("scaleacc_25"))
pred_short <- ggpredict(mo_tu_simp_tf, terms = c("scaleacc_25"))

(graph_tf <- ggplot() +
    geom_line(data = pred_full, color = "red", aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_full, color = "red", aes(ymin = conf.low, ymax = conf.high, 
                                                  x = x), alpha = 0.1, fill = "red") +
    geom_line(data = pred_short, color = "green", aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_short, color = "green", aes(ymin = conf.low, ymax = conf.high, 
                                                     x = x), alpha = 0.1, fill = "green") +
    theme_clean() +
    labs(x = "\nAccessibility", y = "Turnover\n"))
