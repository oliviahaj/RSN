# 1) Cleaning John Yarie's tree ring data to create a file with all of the trees and pith years 
# 2) Creating an updated file with all of the CAFI tree rings available
# 3) Creating a summarized file of CAFI plot ages - across all of 
# - this file will correspond to ages identified in Melissa Boyd's file, 
#   so there will be some trees that are removed from calculating a site average
# - this will also include the fire data where relevant from Adrianna

## OLH
## June 30

# load libraries
library(tidyverse)
library(readxl)

# 1) read in data
## Subplot information on the CAFI sites; this is used to create plot numbers in order
subs <- read.csv("/Users/olhajek/Desktop/RSN/CAFI_data/452_CAFI_SUBPLOT_v20241220.csv")

## This is the age data compiled from Melissa and can be used as a check
mel_age <- read_excel("/Users/olhajek/Desktop/RSN/CAFI_data/CAFI_Boyd_061522/CAFI_stand_age/age_cafi_021621.xlsx", sheet = 1)

## This file contains the tree core data and a page summarizing the ages for each tree, using the summary tab 
trees <- read_excel('/Users/olhajek/Desktop/RSN/CAFI_data/CAFI_Boyd_061522/CAFI_stand_age/AK_PSP_increment_cores USU Mar 2013b.xlsx', sheet = 1, skip = 1) %>%
  filter(if_any(everything(), ~ !is.na(.)))

# 2) Clean up files so that they can be merged
## Extracting a single subplot number for each of the tree; there are a bunch of different formats
tree.2 <- trees %>%
  mutate(sub1 = substr(PSP, 1,3),
    subp = as.numeric(str_extract(sub1, "^\\d+")))

## Creating a matching column between teh PSP on Yarie's file and the subplot CAFI file Miho updated on the BNZ repository
unique.subs <- subs %>%
  select(PLOT, SUBP)%>%
  distinct() %>%
  mutate(subp = row_number())

## For matching purposes, updating the species names to be consistent iwth thsoe in Yarie's tree ring file
mel.age <- mel_age %>%
  rename(PLOT = site) %>%
  # mutate(Spp = case_when(
  #   Sp_age_measured == "birch"  ~ "benteo",
  #   Sp_age_measured == "white spruce"  ~ "picgla",
  #   Sp_age_measured == "aspen"  ~ "poptre",
  #   Sp_age_measured ==  "black spruce"  ~ "picmar",
  #   Sp_age_measured ==  "poplar"  ~ "popbal",
  #   TRUE ~ 'NA'
  # )) %>%
  select(-Sp_age_measured) %>%
  select(c(PLOT, CAFI_sampling_time))


# 3) Join Yaries data with the plot information and create a file that can be merged with the tree_age_compiled
cafi.cores.yarie <- tree.2 %>%
  left_join(unique.subs) %>%
  # select plot, cycle, eyar sampled, subplot, sample number, species, TreeYear, AGE, PithPresent, Notes, Source, Windendro_Compiler
  select(PLOT, SUBP, Year, Spp, `Core Label`, `Count/Years`, Notes) %>%
  rename(YearSampled = Year, SubPlot = SUBP, SampleNumber = `Core Label`, Species = Spp, AGE = `Count/Years`) %>%
  mutate(Species = case_when(
    Species == "benteo" ~ "birch",
    Species == "picgla" ~ "white spruce" ,
    Species == "poptre"~ "aspen",
    Species ==  "picmar"~ "black spruce",
    Species ==  "popbal" ~ "poplar",
    TRUE ~ 'NA'
  )) %>%
  mutate(YearSampled = YearSampled + 2000, TreeYear = YearSampled - AGE + 1, PithPresent = NA, 
         Source = "AK_PSP_increment_cores USU Mar 2013b" ,Windendro_Compiler = "John Yarie, Melissa Boyd",
         # update age to be for 2026?
         AGE = 2026 - TreeYear)  %>%
  # add notes to the trees not included in teh averaging
# Need to join Melissa data to get the Cycle - get rid of the T
  left_join(mel.age) %>%
  mutate(Cycle = substr(CAFI_sampling_time, 2,2)) %>%
  select(-CAFI_sampling_time)%>%
  # ID the trees that will be remove in analysis
  mutate(Notes = case_when(
    SampleNumber == "166-094" ~ "Removed in age calc",
    SampleNumber == "305-094" ~ "Removed in age calc",
    SampleNumber == "088-095" ~ "Removed in age calc",
    SampleNumber == "087-095" ~ "Removed in age calc",
    SampleNumber == "294-370" ~ "Removed in age calc",
    SampleNumber == "316-370" ~ "Removed in age calc",
    SampleNumber == "327-370" ~ "Removed in age calc",
    TRUE ~ Notes
  )) %>%
  # included or removed
  mutate(Remove = ifelse(SampleNumber %in% c("166-094", "305-094","088-095","087-095", "294-370","316-370", "327-370"), 1, 0)) %>%
  rename(Plot = PLOT)

# export this file to be joined
write.csv(cafi.cores.yarie, "/Users/olhajek/Desktop/RSN/CAFI_data/CAFI_Boyd_061522/CAFI_stand_age/cafi_cores_yarie_aged.csv", row.names = FALSE)


# 4) Join with Zach's merged CAFI file that also has Isabelle and Andrew's Cores

cafi.cores <- read.csvcafi.cores <- read.csv("/Users/olhajek/Desktop/RSN/Tree_Age_Data/Compiled/CAFI_Tree_Age_Compiled.csv")

# 5) Update the CAFI file to remove the current Melissa data and replace it with the new tree cores from John Yarie
cafi.cores2 <- cafi.cores %>%
  filter(Source != "age_cafi_021621") %>%
  mutate(Remove = 0) %>%
  rbind(cafi.cores.yarie) %>%
  # Give Andrew's a tree year because they were all sampled in 2025, but going to make it 2026 for consistency with Boyd
  mutate(TreeYear = ifelse(Source == "CAFI_Tree_Cores_2025", 2025 - AGE +1, TreeYear))

# 6) Save this compiled file with all of teh tree ages
write.csv(cafi.cores2, "/Users/olhajek/Desktop/RSN/Tree_Age_Data/Compiled/CAFI_Tree_Age_Compiled_OH.csv", row.names = FALSE)

# 7) Use this file to create an update, compiled CAFI age data
# - want to compare ages across sampling cores/efforts and see what sort of overlap there is
str(cafi.cores2)

cafi.ages <- cafi.cores2 %>%
  mutate(Remove = ifelse(Notes == "Strongly recommend omitting from data. Very fractured, P2 crosses through bark and extensive scar, not counted (or final segments connected)",
                         1, Remove)) %>%
  mutate(Remove = ifelse(is.na(Remove), 0, 1))%>%
  filter(Remove != 1) %>%
  group_by(Plot, YearSampled, Windendro_Compiler) %>%
  summarize(n = n(), mean_ty = mean(TreeYear, na.rm = TRUE), min_ty = min(TreeYear), max_ty = max(TreeYear), 
            mean_AGE = mean(AGE, na.rm = TRUE), median_AGE = median(AGE, na.rm = TRUE), min_AGE = min(AGE), max_AGE = max(AGE)) %>%
  mutate(range_AGE = max_AGE - min_AGE, mean_med = mean_AGE - median_AGE) %>%
  filter(!is.na(Plot))

# How many plots have multiple sampling efforts
sample.num <- cafi.ages %>%
  ungroup() %>%
  select(c(Plot, Windendro_Compiler)) %>%
  unique() %>%
  group_by(Plot) %>%
  summarize(n=n())
# only 2 plots were sampled multiple times - 1022 and 1088
# 1088 is quite different - see if it burned?
# check with the raw data on these...see if there is any indication of fire? can look at trees


# 8) Save the plot ages summary
write.csv(cafi.ages, "/Users/olhajek/Desktop/RSN/Tree_Age_Data/Compiled/CAFI_Ages_Summarized.csv",row.names = FALSE)
