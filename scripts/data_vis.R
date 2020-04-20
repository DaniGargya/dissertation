# Data visualisation
# Dani Gargya
# March 2020

# Load data ----
biotime_full <- read.csv("data/BioTIMEQuery02_04_2018.csv")
biotime_meta <- read.csv("data/BioTIMEMetadata_02_04_2018.csv")
bio <- read.csv("data/bio.csv")

bio_short <- bio %>% 
  distinct(STUDY_ID_PLOT, STUDY_ID, CENT_LAT, CENT_LONG, TAXA, TOTAL, START_YEAR, END_YEAR, NUMBER_OF_SAMPLES, AREA_SQ_KM, HAS_PLOT, HAS_PLOT, NUMBER_LAT_LONG)

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

(graph_acc <- ggplot() +
  geom_line(data = predictions_acc, aes(x = x, y = predicted),
            size = 2) +
  geom_ribbon(data = predictions_acc, aes(ymin = conf.low, ymax = conf.high, 
                                        x = x), alpha = 0.1) +
  geom_point(data = data1, aes(x = scaleacc_25, y = Jtu),
             alpha = 0.1, size = 2) +
  #annotate("text", x = -0.65, y = 5, label = "Slope = -0.06, Std. error = 0.01") +  
  #scale_x_continuous(limits = c (0.8, 1)) +
  theme_clean() +
  labs(x = "\nAccessibility", y = "Jaccard dissimilarity\n"))

ggsave(graph_acc, filename = "outputs/graph_acc.png",
       height = 5, width = 8)


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





#### RQ2: taxa making raincloud plot ----
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
me <- ggpredict(mo_tu2, terms = c("scaleacc_25", "TAXA"), type = "re")

# plot
(graph_acc <- ggplot() +
   geom_line(data = me, aes(x = x, y = predicted),
             size = 2) +
   geom_ribbon(data = me, aes(ymin = conf.low, ymax = conf.high, 
                                           x = x), alpha = 0.1) +
   facet_wrap("group") +
   geom_point(data = data1, aes(x = scaleacc_25, y = Jtu),
              alpha = 0.1, size = 2) +
   facet_wrap ("TAXA") +
   #annotate("text", x = -0.65, y = 5, label = "Slope = -0.06, Std. error = 0.01") +  
   #scale_x_continuous(limits = c (0.8, 1)) +
   theme_clean() +
   labs(x = "\nAccessibility", y = "Jaccard dissimilarity\n"))

plot(me)

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

# from g ----

pd <- position_dodge(0.2) # So that the error bars on graphs don't overlap 

(pred_plot <- ggplot(predictSP, 
                    aes(x=TAXA, y=mean, colour=treatment, group=treatment))+ 
  geom_errorbar(aes(ymin=down, ymax=up), 
                colour="black", width=.2, position=pd) + 
  geom_point(position=pd, size=4) + 
  theme_classic() + 
  labs(x="Taxa", y="Jaccard dissimilarity") + 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14)))


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
  geom_line(data = predictions_2, aes(x = x, y = predicted),
            size = 2) +
  geom_ribbon(data = predictions_2, aes(ymin = conf.low, ymax = conf.high, 
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
pre_low <- predictions_3 %>% 
  filter(hpd == "Low")

pre_mod <- predictions_3 %>% 
  filter(hpd == "Moderate")

pre_high <- predictions_3 %>% 
  filter(hpd == "High")



# filter dataset all ----
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


load("outputs/IMSsimple_model1.RData")
summary(mo_tu_simp1)

