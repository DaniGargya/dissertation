# Sensitivity analysis dissertation
# Dani Gargya
# April 2020

update.packages(ask = FALSE, checkBuilt = TRUE)  # update R packages
tinytex::tlmgr_update() 
# to run
# scale 100 with 2 random effects
# smaller time frame
# richness trends
# area?

# data ----
bio <- read.csv("data/bio.csv") %>%  dplyr::select(-X)

# sites coincidence with protected areas ----
protected_area <- data1_plants %>% 
  group_by(PROTECTED_AREA) %>% 
  summarise(plots =length(unique(STUDY_ID_PLOT)))
# 32.6% fall into protected areas
# plants: 38.9%
# birds: 0.01%
# terrestrial invertebrates: 44.5%
# mammals: 19.8

# smaller timeframe ----
short_timeframe <- bio %>% 
  filter(START_YEAR > 1980,
         END_YEAR < 2015)

# scale accessibility ----
data1 <- data1 %>% 
  mutate(scaleacc_100= 1 - ((e_100-min(e_100))/(max(e_100)-min(e_100)))) %>% 
  mutate(scaleacc_50= 1 - ((e_50-min(e_50))/(max(e_50)-min(e_50)))) %>% 
  mutate(scalehpd_100= (e_hpd100-min(e_hpd100))/(max(e_hpd100)-min(e_hpd100))) %>% 
  mutate(scalehpd_50= (e_hpd50-min(e_hpd50))/(max(e_hpd50)-min(e_hpd50))) %>% 
  left_join(data1_1, by = "STUDY_ID_PLOT")

check <- data1 %>% 
  drop_na(e_hpd)
min(check$e_hpd)
max(check$e_hpd)

hist(data1$scalehpd_100)
hist(data1$scaleacc_25)


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
summary(mo_tu_scale100)

# sensitivity analysis scale (50km)
mo_tu_simp7 <- brm(bf(Jtu ~ scaleacc_50 + scalehpd_50 + duration_plot + TAXA +
                        (1|cell) + (1|STUDY_ID)), # or with (1|cell)?
                   family = zero_one_inflated_beta(), 
                   data = data1,
                   iter = 4000,
                   warmup = 1000,
                   inits = '0',
                   control = list(adapt_delta = 0.85),
                   cores = 2, chains = 4)

save(mo_tu_simp7, file = "outputs/mo_tu_simp7.RData")

# differences in latitude? (more promminent in tropics?)----

# number of data points?

# model only terrestrial plants ----
# sensitivity analysis plants
data1_plants <- data1 %>% 
  filter(TAXA == "Terrestrial plants")

mo_tu_simp_plants <- brm(bf(Jtu ~ scaleacc_25 + scalehpd_25 + duration_plot +
                              (1|cell) + (1|STUDY_ID)), 
                         family = zero_one_inflated_beta(), 
                         data = data1_plants,
                         iter = 4000,
                         warmup = 1000,
                         inits = '0',
                         control = list(adapt_delta = 0.85),
                         cores = 2, chains = 4)

save(mo_tu_simp_plants, file = "outputs/mo_tu_simp_plants.RData")

# smaller time frame 1970 - 2010 ----
data1_tfsmall <- data1 %>% 
  filter(!START_YEAR < 1970) %>% 
  filter(!END_YEAR > 2015)

mo_tu_simp_tf <- brm(bf(Jtu ~ scaleacc_25 + scalehpd_25 + duration_plot + TAXA +
                              (1|cell) + (1|STUDY_ID)), 
                         family = zero_one_inflated_beta(), 
                         data = data1_tfsmall,
                         iter = 4000,
                         warmup = 1000,
                         inits = '0',
                         control = list(adapt_delta = 0.85),
                         cores = 2, chains = 4)

save(mo_tu_simp_tf, file = "outputs/mo_tu_simp_tf.RData")
summary(mo_tu_simp_tf)
# richness trends ----
