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

# differences in latitude? (more promminent in tropics?)----

# number of data points?

# model only terrestrial plants ----

