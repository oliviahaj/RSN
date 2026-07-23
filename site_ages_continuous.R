# Continuous age file
## The purpose of this script is to take the inventory data and identify all of the unique sample periods
## for each site and assign a year adn age to that year
## this is because sites change over time. Also includean indication of how many plots were sampled for 
## a given year. this should already be available in the RSN site_years.csv

## OLH
## July 23, 2026

##############################################
# 1) Load libraries
library(tidyverse)

##############################################
# 2) Read in data
# Inventory Data

trees <- read.csv("/Users/olhajek/Desktop/RSN/RSN_proj/Data/harmonized/RSN_CAFI_merged.csv")

# Tree Age file - combined RSN and CAFI
# Note that this first file is for those with tree ages adn pith year! it does not have the sites with a burn year in them

ages <- read.csv("/Users/olhajek/Desktop/RSN/RSN_proj/Data/Tree_Age_Data/Compiled/RSN_CAFI_Summarized_Ages.csv")

# Site summary - burn year for RSN sites
burns <- read.csv('/Users/olhajek/Desktop/RSN/RSN_proj/Data/2026 Data/RSN_Site/RSN_Master_Site.csv')

# Site years
# RSN
site_yrs_rsn <- read.csv('/Users/olhajek/Desktop/RSN/RSN_proj/Data/2026 Data/Data_checks/site_years.csv')

# CAFI
site_yrs_cafi <- read.csv('/Users/olhajek/Desktop/RSN/RSN_proj/Data/2026 Data/CAFI/452_CAFI_SUBPLOT_v20241220.csv')

##############################################
# 3) Process data files and merge everything
# For the inventory data, select each unique site-year-number of subplots

glimpse(trees)
trees.summary <- trees %>%
  group_by(Year, PLOT, Network) %>%
  summarize(
    # no of unique subplots
    no_subp = length(unique(SUBP)))
  
# prep the site-yrs to be joined
glimpse(site_yrs_rsn)
rsn <- site_yrs_rsn %>%
  select(-c(X, age_class_yr)) %>%
  rename(PLOT = SITE)

glimpse(site_yrs_cafi)
cafi <- site_yrs_cafi %>%
  select(PLOT, CYCLE, SUBP, MEASYEAR, NOTES)%>%
  rename(Year = MEASYEAR,
         sample.rd = CYCLE) %>%
  group_by(PLOT, Year, sample.rd) %>%
  summarize(n_plot = length(unique(SUBP)),
            notes_samplingarea = if (all(is.na(NOTES))) NA_character_ else paste(NOTES[!is.na(NOTES)], collapse = "; ")) %>%
  ungroup()%>%
  mutate(sample.yr = NA, 
         sampled_plots = NA,
         sampled_plots_period = n_plot)

site_yrs <- rbind(rsn, cafi)

# Merge the site-yrs with the 



##############################################
# 2) Export data file





