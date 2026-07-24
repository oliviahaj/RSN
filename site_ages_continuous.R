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

# Recent burns
rb <- read.csv("/Users/olhajek/Desktop/RSN/RSN_proj/Data/Tree_Age_Data/Source/RSN/recent_RSN_burns.csv")

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

# Merge the burns with the ages
glimpse(burns)
glimpse(ages)

burns.2 <- burns %>%
  select(sitename, burn_year) %>%
  filter(!is.na(burn_year))%>%
  rename(PLOT = sitename) %>%
  mutate(burn_AGE = 2026-burn_year) 

ages.2 <- ages %>%
  rename(PLOT = Plot) %>%
  select(-c(median_AGE, min_AGE, max_AGE, range_AGE, mean_med))

age.combo <- full_join(burns.2, ages.2) %>%
  # check difference between mean age and the burn year
# max difference is 22 years with some of the gerstel river sites; nwill have both, but 
# don't need to have teh difference in the file to combine
  mutate(diff = burn_year - mean_ty)
  

# Merge the site-yrs with the age file
glimpse(site_yrs)
glimpse(age.combo)

site_yr_ages <- left_join(site_yrs, age.combo)

# Update this with the recent burns file
glimpse(rb)
glimpse(site_yr_ages)

# Merge this with the tree files - which plots and years, etc are missing data

sya <- site_yr_ages %>%
  left_join(rb) %>%
  # going to get rid of the burn year for FP5C for just the recent years
  mutate(burn_year = ifelse(PLOT == "FP5C", NA, burn_year)) %>%
  # going to update the burn year when it is greater
  mutate(burn_year = ifelse(!is.na(burn_year2) & Year > burn_year2, burn_year2, burn_year), 
         future_year = ifelse(Year < burn_year2,  burn_year2, NA)) %>%
  # get rid of unnecessary age columsn 
  select(-c(mean_AGE, diff, burn_AGE, burn_year2)) %>%
  # going to select an age based on the difference between the years
  mutate(year_diff = burn_year - floor(mean_ty)) %>%
  # going to coalesce on this to select an age, if greater than 10 years diff based on above
  # then will select an age
  mutate(age1 = ifelse(year_diff>10, burn_year, mean_ty), 
         age = floor(coalesce(age1, mean_ty, burn_year))) %>%
  select(-c(age1, year_diff))


glimpse(trees.summary)
glimpse(sya)

# creating a new row for the ages because it is missing in the CAFI subplot file
sya.1011 <- sya %>%
  filter(PLOT == "1011") %>%
  filter(Year == 2015) %>%
  mutate(Year = 2020)

# add it back in
sya <- rbind(sya, sya.1011)


# Join the site-year ages with the tree summary
tree.ages <- left_join(trees.summary, sya) %>%
  # chekcing to see which plots are aligned
 # mutate(sample.diff = n_plot - no_subp, sample.diff2 = sampled_plots - no_subp)
  # going to remove some of hte different subplot number and sampling things
  select(-c(n_plot))
  
# check the NA ages - which plots don't have any age
na.ages <- tree.ages %>%
  filter(is.na(age))


cafi.sites <- trees.summary %>%
  ungroup()%>%
  filter(Network == "CAFI") %>%
  select(PLOT)%>%
  distinct()


# create a new sample.pd column? 
unique(tree.ages$Year)

ggplot(tree.ages, aes(Year, PLOT))+
  geom_point()+
  facet_wrap(~Network, nrow=1, scales="free")

tree.ages <- tree.ages %>%
  mutate(sample.pd = case_when(
    Year < 2016 ~ 1, 
    Year > 2019 ~ 3, 
    TRUE ~ 2
  ))

##############################################
# 2) Export data file
# I think we want the file with teh data pre-tree, so site_yr_ages with everythign that we have

write.csv(tree.ages, "/Users/olhajek/Desktop/RSN/RSN_proj/Data/harmonized/continuous_tree_ages.csv")






