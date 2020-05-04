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

# smaller time frame
load("outputs/mo_tu_simp_tf.RData")
summary(mo_tu_simp_tf)

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
display.brewer.pal(n=4, name = "Set1")
brewer.pal(n=4, name = "Set1")
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



# panel time-series ----

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

(acc_d <- ggMarginal(graph_acc, type="density", size = 3, fill = "#666666", xparams = list(fill = "#A6761D")))

ggsave(acc_d, file = "outputs/graph_acc_den.png" , width = 7, height = 5)

#### RQ2: turnover ~ hpd ----
# prediction
pred_hpd <- ggpredict(mo_tu_simp6, terms = c("scalehpd_25"))

# model vis
(graph_hpd <- ggplot() +
  geom_line(data = pred_hpd, aes(x = x, y = predicted),
            size = 2, color = "#A6761D") +
  geom_ribbon(data = pred_hpd, aes(ymin = conf.low, ymax = conf.high, 
                                     x = x), fill = "#A6761D", alpha = 0.1) +
  geom_point(data = data1, aes(x = scalehpd_25, y = Jtu),
             alpha = 0.1, size = 2) +
  theme_clean() +
  labs(x = "\nHuman population density (proportion)", y = "Turnover\n") +
  theme(legend.position = c(0.85, 0.8)))

ggsave(graph_hpd, file = "outputs/graph_hpd.png", width = 7, height = 5)

(hpd_d <- ggMarginal(graph_hpd, type="density", size = 3, fill = "#666666", xparams = list(fill = "#A6761D")))

ggsave(hpd_d, file = "outputs/graph_hpd_den.png" , width = 7, height = 5)

#### RQ3: turnover ~ taxa ----
# raincloud plot
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
    geom_point(aes(size = 1)) +
    geom_pointrange(aes(ymin = conf.low, ymax =  conf.high), size = 2) +
    scale_fill_manual(values = taxa.palette) +
    scale_colour_manual(values = taxa.palette) +
    theme_clean() +
    theme(legend.position = "none",
          axis.text.y=element_blank(),
          axis.title.y=element_blank()) +
    geom_hline(yintercept=0, linetype="dashed", size=1) +
    labs(x = NULL, y = "Effect size\n") +
    coord_flip())

ggsave(graph_taxa, filename = "outputs/taxa_model.png",
       height = 5, width = 5)

panel_taxa <- ggarrange(raincloud_taxa, graph_taxa, labels = c("A", "B"), ncol = 2)
panel_taxa2 <- grid.arrange(raincloud_taxa, graph_taxa, respect=TRUE, ncol = 2,  widths=c(1.5,1))

#graph_taxa = graph_taxa + theme(aspect.ratio=2)

#p <- grid.arrange(raincloud_taxa, graph_taxa,  ncol = 2) # one is square, the other thinner  
  

ggsave(p, filename ="outputs/panel_taxa3.png",
       height = 5, width = 8)

#### sensitivtiy analysis ----
# plants vs all ----
pred_plants <- ggpredict(mo_tu_simp_plants, terms = c("scaleacc_25"))
pred_all <- ggpredict(mo_tu_simp6, terms = c("scaleacc_25"))

data1_plants <- data1 %>% 
  filter(TAXA == "Terrestrial plants")

data1_nopl <- data1 %>% 
  filter(!TAXA == "Terrestrial plants")

(graph_plants <- ggplot() +
    geom_line(data = pred_plants, color = "#E6AB02", aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_plants, color = "#E6AB02", aes(ymin = conf.low, ymax = conf.high, 
                                                          x = x), alpha = 0.1, fill = "#E6AB02") +
    geom_point(data = data1_plants, color = "#E6AB02", aes(x = scaleacc_25, y = Jtu),
               alpha = 0.1, size = 2) +
    geom_line(data = pred_all, color = "#66A61E", aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_all, color = "#66A61E", aes(ymin = conf.low, ymax = conf.high, 
                                                        x = x), alpha = 0.1, fill = "#66A61E") +
    geom_point(data = data1_nopl, color = "#66A61E", aes(x = scaleacc_25, y = Jtu),
               alpha = 0.1, size = 2) +
    theme_clean() +
    labs(x = "\nAccessibility", y = "Turnover\n"))

ggsave(graph_plants, filename = "outputs/graph_plants.png",  height = 5, width = 8)

# scales ----
pred_1 <- ggpredict(mo_tu_simp8, terms = c("scaleacc_1"))
pred_25 <- ggpredict(mo_tu_simp6, terms = c("scaleacc_25"))
pred_50 <- ggpredict(mo_tu_simp7, terms = c("scaleacc_50"))
pred_100 <- ggpredict(mo_tu_scale100, terms = c("scaleacc_100"))

(graph_scale <- ggplot() +
    geom_line(data = pred_1, color = "#E41A1C", aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_1, color = "#E41A1C", aes(ymin = conf.low, ymax = conf.high, 
                                                  x = x), alpha = 0.1, fill = "#E41A1C") +
    geom_line(data = pred_25, color = "#377EB8" , aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_25, color = "#377EB8" , aes(ymin = conf.low, ymax = conf.high, 
                                                     x = x), alpha = 0.1, fill = "#377EB8" ) +
    geom_line(data = pred_50, color = "#4DAF4A" , aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_50, color = "#4DAF4A", aes(ymin = conf.low, ymax = conf.high, 
                                                    x = x), alpha = 0.1, fill = "#4DAF4A" ) +
    geom_line(data = pred_100, color = "#984EA3", aes(x = x, y = predicted),
              size = 2) +
    geom_ribbon(data = pred_100, color = "#984EA3", aes(ymin = conf.low, ymax = conf.high, 
                                                       x = x), alpha = 0.1, fill = "#984EA3") +
    theme_clean() +
    labs(x = "\nAccessibility", y = "Turnover\n") +
    ylim(0, 1))

ggsave(graph_scale, filename = "outputs/graph_scale.png", height = 5, width = 8)
#### other extra graphs ----
# raincloud plot taxa accessibility ----
(raincloud_acc <- 
   ggplot(data = data1, 
          aes(x = reorder(TAXA, desc(scaleacc_25)), y = scaleacc_25, fill = TAXA)) +
   geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.8) +
   geom_point(aes(y = scaleacc_25, color = TAXA), 
              position = position_jitter(width = 0.15), size = 1, alpha = 0.1) +
   geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
   labs(y = "\nAccessibility (proportion)", x = NULL) +
   guides(fill = FALSE, color = FALSE) +
   scale_fill_manual(values = taxa.palette) +
   scale_colour_manual(values = taxa.palette) +
   coord_flip() +
   theme_clean())

ggsave(raincloud_acc, filename = "outputs/raincloud_taxa_acc_graph.png",
       height = 5, width = 8)

# protected areas ----
(graph_proteccted_areas <- ggplot() +
   geom_point(data = data1, aes(x = scaleacc_25, y = Jtu, colour = PROTECTED_AREA),
              alpha = 0.8, size = 1.5) +
   scale_colour_manual(values = c("#E6AB02", "#666666")) +
   theme_clean() +
   theme(legend.position = "none") +
   labs(x = "\nAccessibility", y = "Turnover\n"))

ggsave(graph_proteccted_areas, filename = "outputs/graph_protected_areas.png",
       height = 5, width = 8)

# graph duration and turnover ----
pred_duration <- ggpredict(mo_tu_simp6, terms = c("duration_plot"))

(graph_duration <- ggplot() +
    geom_line(data = pred_duration, aes(x = x, y = predicted),
              size = 2, color = "#666666") +
    geom_ribbon(data = pred_duration, aes(ymin = conf.low, ymax = conf.high, 
                                     x = x), fill = "#666666", alpha = 0.1) +
    theme_clean() +
    labs(x = "\nMonitoring duration (in years)", y = "Turnover\n") +
    ylim(0, 1))

ggsave(graph_duration, file = "outputs/graph_duration.png", width = 7, height = 5)

# histogram distribution of random and real values ----
# accessibility
bio_short <- bio %>% 
  distinct(STUDY_ID_PLOT, .keep_all = TRUE)

SP <- bio_short %>% 
  dplyr::select(LATITUDE, LONGITUDE, STUDY_ID_PLOT) %>% 
  distinct(LATITUDE, .keep_all = TRUE)

library(generator)
fake_lat <- r_latitudes(1023)
fake_long <- r_longitudes(1023)

fake_ll <- SP %>% 
  mutate(fake_lat = c(fake_lat),
         fake_long = c(fake_long)) %>% 
  dplyr::select(- LATITUDE, -LONGITUDE)

f_points <- cbind(fake_ll$fake_long, fake_ll$fake_lat)


f_sppoints <- SpatialPoints(f_points, proj4string=CRS('+proj=longlat +datum=WGS84'))
f_tp <- spTransform(f_sppoints, crs(aa))

f_e <- extract(aa, f_tp)

f_bio_aa <- cbind(fake_ll, f_e)
f_bio_aa_short <- na.omit(f_bio_aa)

# hpd
f_tp_hpd <- spTransform(f_sppoints, crs(hpd))

f_e_hpd <- extract(hpd, f_tp_hpd)

f_bio_hpd <- cbind(fake_ll, f_e_hpd)
f_bio_hpd_short <- na.omit(f_bio_hpd)

# figure overall
png("outputs/fake_ll.png")
par(
  mfrow=c(1,2),
  mar=c(4,4,1,0)
)
hist(f_bio_aa_short$f_e, col=rgb(1,0,0,0.5), xlab="Accessibility score", ylab ="Number of points", main = "") # more normal distributed?
hist(data1$e, col=rgb(0,0,1,0.5), add=T)
hist(f_bio_hpd_short$f_e_hpd, col=rgb(1,0,0,0.5), xlab="Human population density score", ylab ="Number of points", main = "")
hist(data1$e_hpd, col=rgb(0,0,1,0.5), add=T)
legend("topright", legend=c("Random","BioTIME"), col=c(rgb(1,0,0,0.5), 
                                                       rgb(0,0,1,0.5)), pt.cex=2, pch=15)
dev.off()

# histogram acc and hpd ----
png("outputs/hist_acc_hpd.png")
hist(data1$scaleacc_25, col = "#1B9E77", xlab="Accessibility/Human population density", ylab ="Number of points", main = "")
hist(data1$scalehpd_25, col = "#A6761D", add =T)
legend("top", legend=c("Accessibility","Human population density"), col=c("#1B9E77", "#A6761D"), pt.cex=2, pch=15)
dev.off()      
