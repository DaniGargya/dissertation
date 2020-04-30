# Data vis clean
# Dani Gargya
# April 2020

# Load data ----
data1 <- read.csv("data/data1.csv") %>%  dplyr::select(-X)

# MODELS TO USE
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

# create colour palette ----
display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n = 8, name = "Dark2")
taxa.palette <- c("#D95F02", "#7570B3", "#E7298A", "#E6AB02")
names(taxa.palette) <- levels(data1$TAXA)
# colours turnover c("#7CFC00")
c("#CD9B1D", "#FFC125", "#FFFFFF")

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



#### RQ1: turnover ~ accessibility + duration ----
# prediction
pred_acc_d <- ggpredict(mo_tu_simp6, terms = c("scaleacc_25", "duration_plot"))

pred_acc_d$Monitoring <- factor(pred_acc_d$group, levels = c("6", "19.1", "32.2"),
                                 labels = c("Short-term", "Moderate", "Long-term"))

# calculate change
pred_acc <- ggpredict(mo_tu_simp6, terms = c("scaleacc_25"))

# model vis
(graph_acc <- ggplot() +
  geom_line(data = pred_acc_d, aes(x = x, y = predicted, color = Monitoring),
            size = 2) +
  geom_ribbon(data = pred_acc_d, aes(ymin = conf.low, ymax = conf.high, 
                                         x = x, fill = Monitoring), alpha = 0.1) +
  geom_point(data = data1, aes(x = scaleacc_25, y = Jtu),
             alpha = 0.1, size = 2) +
  theme_clean() +
  scale_fill_manual(values = c("#7CFC00", "#66A61E", "#1B9E77")) +
  scale_colour_manual(values = c("#7CFC00", "#66A61E", "#1B9E77")) +
  labs(x = "\nAccessibility (proportion)", y = "Turnover\n") +
  theme(legend.position = c(0.2, 0.8)))

ggsave(graph_acc, file = "outputs/graph_acc.png", width = 7, height = 5)

(acc_d <- ggMarginal(graph_acc, type="density", size = 3, fill = "green"))

ggsave(acc_d, file = "outputs/graph_acc_den.png" , width = 7, height = 5)

#### RQ2: turnover ~ hpd + duration ----
# prediction
pred_hpd_d <- ggpredict(mo_tu_simp6, terms = c("scalehpd_25", "duration_plot"))

pred_hpd_d$Monitoring <- factor(pred_hpd_d$group, levels = c("6", "19.1", "32.2"),
                                labels = c("Short-term", "Moderate", "Long-term"))

# model vis
(graph_hpd <- ggplot() +
  geom_line(data = pred_hpd_d, aes(x = x, y = predicted, color = Monitoring),
            size = 2) +
  geom_ribbon(data = pred_hpd_d, aes(ymin = conf.low, ymax = conf.high, 
                                     x = x, fill = Monitoring), alpha = 0.1) +
  geom_point(data = data1, aes(x = scalehpd_25, y = Jtu),
             alpha = 0.1, size = 2) +
  theme_clean() +
  scale_fill_manual(values = c("#FFC125" , "#CD9B1D", "#A6761D")) +
  scale_colour_manual(values = c("#FFC125" , "#CD9B1D", "#A6761D")) +
  labs(x = "\nHuman population density (proportion)", y = "Turnover\n") +
  theme(legend.position = c(0.85, 0.8)))

ggsave(graph_hpd, file = "outputs/graph_hpd.png", width = 7, height = 5)

(hpd_d <- ggMarginal(graph_hpd, type="density", size = 3, fill = "brown"))

ggsave(hpd_d, file = "outputs/graph_hpd_den.png" , width = 7, height = 5)

#### RQ3: turnover ~ taxa ----
# raincloud plot ----
# by Ben Marwick
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

# raincloud plot taxa jtu
(raincloud_taxa <- 
   ggplot(data = data1, 
          aes(x = reorder(TAXA, desc(Jtu)), y = Jtu, fill = TAXA)) +
   geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.8) +
   geom_point(aes(y = Jtu, color = TAXA), 
              position = position_jitter(width = 0.15), size = 1, alpha = 0.1) +
   geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
   labs(y = "\n Turnover", x = NULL) +
   guides(fill = FALSE, color = FALSE) +
   scale_y_continuous(limits = c(0, 1)) +
   scale_fill_manual(values = taxa.palette) +
   scale_colour_manual(values = taxa.palette) +
   coord_flip() +
   theme_clean())

ggsave(raincloud_taxa, filename = "outputs/raincloud_taxa_graph.png",
       height = 5, width = 8)


# effect size graph
pred_taxa <- ggpredict(mo_tu_simp6, terms = c("TAXA"))

(graph_taxa <- ggplot(pred_taxa, aes(x = reorder(x, desc(predicted)), y = predicted, color = x)) +
    geom_point(aes(size = 2)) +
    geom_pointrange(aes(ymin = conf.low, ymax =  conf.high), size = 2) +
    scale_fill_manual(values = taxa.palette) +
    scale_colour_manual(values = taxa.palette) +
    theme_clean() +
    theme(legend.position = "none",
          axis.text.x=element_blank(),
          axis.title.x=element_blank()) +
    geom_hline(yintercept=0, linetype="dashed", size=1) +
    labs(x = NULL, y = "Effect size\n") +
    coord_flip())

panel_taxa <- ggarrange(raincloud_taxa, graph_taxa, ncol = 2)

ggsave(panel_taxa, filename ="outputs/panel_taxa.png",
       height = 5, width = 8)
