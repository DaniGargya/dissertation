# Sensitivity analysis dissertation
# Dani Gargya
# April 2020

# data ----
bio <- read.csv("data/bio.csv") %>%  dplyr::select(-X)

# sites coincidence with protected areas ----
protected_area <- bio %>% 
  group_by(PROTECTED_AREA) %>% 
  summarise(plots =length(unique(STUDY_ID_PLOT)))

# smaller timeframe ----
short_timeframe <- bio %>% 
  filter(START_YEAR > 1980,
         END_YEAR < 2015)

# scale accessibility ----
data1 <- data1 %>% 
  mutate(scaleacc_100= 1 - ((e_100-min(e_100))/(max(e_100)-min(e_100)))) %>% 
  mutate(scalehpd_100= (e_hpd100-min(e_hpd100))/(max(e_hpd100)-min(e_hpd100)))

check <- data1 %>% 
  drop_na(e_hpd)
min(check$e_hpd)
max(check$e_hpd)

hist(data1$scalehpd_100)
hist(data1$scaleacc_25)


mo_tu_scale100 <- brm(bf(Jtu ~ scaleacc_100 + scalehpd_100 + duration_plot + TAXA +
                        (1|STUDY_ID)),
                   family = zero_one_inflated_beta(), 
                   data = data1,
                   iter = 4000,
                   warmup = 1000,
                   inits = '0',
                   control = list(adapt_delta = 0.85),
                   cores = 2, chains = 4)

# differences in latitude? (more promminent in tropics?)----

# number of data points?

# model only terrestrial plants ----

