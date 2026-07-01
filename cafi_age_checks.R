# 1) Cleaning John Yarie's tree ring data to create a file with all of the trees and pith years 
# 2) Creating a summary file of plot ages
# - this file will correspond to ages identified in Melissa Boyd's file, so there will be some trees that are removed from calculating a site average

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
         source = "AK_PSP_increment_cores USU Mar 2013b" ,Windendro_Compiler = "John Yarie, Melissa Boyd",
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
  mutate(Remove = ifelse(SampleNumber %in% c("166-094", "305-094","088-095","087-095", "294-370","316-370", "327-370"), 1, 0))

# export this file to be joined
write.csv(cafi.cores.yarie, "/Users/olhajek/Desktop/RSN/CAFI_data/CAFI_Boyd_061522/CAFI_stand_age/cafi_cores_yarie_aged.csv", row.names = FALSE)


# 4) Filter plots where the ages are not just averages of all trees 
tree.32 <- tree.2 %>%
  left_join(unique.subs)%>%
  filter(PLOT %in% c(1180, 1096, 1175,1183))
# all but 1180 are just an average of the two oldest trees adn then the timepoint. not sure what's going on wiht 1180

tree.3 <- tree.2 %>%
  left_join(unique.subs) %>%
  group_by(PLOT, Spp) %>%
  summarize(mean = mean(`Count/Years`), med = median(`Count/Years`),max = max(`Count/Years`),min = min(`Count/Years`), n = n())%>%
  left_join(mel.age) %>%
  mutate(diff = age_biomass_yr_msrd - mean) %>%
  filter(diff > 1 | diff < -1)

ggplot(tree.3, aes(diff))+
  geom_histogram()+
  facet_wrap(~CAFI_sampling_time)

tree.4 <- tree.2 %>%
  left_join(unique.subs) %>%
  group_by(PLOT) %>%
  summarize(mean = mean(`Count/Years`), max = max(`Count/Years`),min = min(`Count/Years`))%>%
  left_join(mel.age)%>%
  mutate(diff = age_biomass_yr_msrd - mean)


subp <- read.csv("/Users/olhajek/Desktop/RSN/CAFI_data/452_CAFI_SUBPLOT_v20241220.csv")

str(subp)
subp.prep <- subp %>%
  select(PLOT, CYCLE, MEASYEAR) %>%
  mutate(CAFI_sampling_time = paste("T",CYCLE, sep="")) %>%
  distinct()

test <- left_join(tree.3, subp.prep)  

test.2 <- test %>%
  mutate(year2 = Year + 2000, 
         diff2 =  MEASYEAR - year2, 
         off = diff2 - diff, 
         range = max-min)
