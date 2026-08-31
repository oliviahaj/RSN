# RSN Age checks
# Using the data from Zach to check on the ages have and create a compiled file similar to what was created with CAFI

# load libraries
library(tidyverse)

# Read in data
rsn.trees <- read.csv("/Users/olhajek/Desktop/RSN/RSN_proj/Data/Tree_Age_Data/Compiled/RSN_TreeAge_Compiled.csv")

# Update the dataframe for compatability with the CAFI Dataset
# need to add in AGE, rename, plot, subplot
# A few notes on tree ages
glimpse(rsn.trees)

# which sites have more than one data sources
rsn.sources <- rsn.trees %>%
  select(Site, Source, Data_Sheet) %>%
  distinct()

# this code here summarizes by tree first! now going to do by stand across plots
rsn.trees.2 <- rsn.trees %>%
  # get rid of cbsnbc - charred, burned standing
  filter(Species != "CBSNBC" | is.na(Species))%>%
  rename(SubP = Plot, Plot = Site, YearSampled = Year_Collected, TreeYear = Pith_Year) %>%
  # Going to get keep the Jamie found source for WCM1 and DCM1; for WCM2 just going to merge (113 vs. 180, different but maybe not hugely meaningfully)
  # filter(Plot != "WCM1" | Source != "MAW, CJ") %>%
  # filter(Plot != "DCM1" | Source != "MAW, LP") %>%
  group_by(Plot, SubP, Sample_Number, Species, Source) %>%
  #group_by(Plot, YearSampled, Source) %>%
  summarize(n = n(), mean_tree = mean(TreeYear, na.rm = TRUE), min_tree = min(TreeYear), max_tree = max(TreeYear),  
            mean_rc_tree = mean(Ring_Count, na.rm=TRUE), diff = max_tree - min_tree) %>%
  mutate(SPECIES = case_when(
    Species == "LALA" ~ "LARLAR", 
    Species == "LALA" ~ "LARLAR",
    Species == "PICLGLA" ~ "PICGLA",
    Species == "Picmar" ~ "PICMAR",
    Species == "PIGL" ~ "PICGLA",
    Species == "PIMA" ~ "PICMAR",
    TRUE ~ Species
  ))

rsn.trees.3 <- rsn.trees.2 %>%
  ungroup()%>%
  mutate(AGE = 2026 - mean_tree) %>%
  group_by(Plot) %>%
  summarize(n = n(), mean_ty = mean(mean_tree, na.rm = TRUE), min_ty = min(mean_tree), max_ty = max(mean_tree),  mean_RC = mean(mean_rc_tree, na.rm=TRUE),
            mean_AGE = mean(AGE, na.rm = TRUE), median_AGE = median(AGE, na.rm = TRUE), min_AGE = min(AGE), 
            max_AGE = max(AGE), sd_AGE = sd(AGE, na.rm =TRUE)) %>%
  mutate(Notes = case_when(
    Plot == "WCM1" ~ "Kept Jamie found data; seemed to repeat data",
    Plot == "WCM2" ~ "Kept Jamie found data; seemed to repeat data",
    Plot == "DCM1" ~ "Kept Jamie found data; huge difference from other data file; confirm",
    Plot == "FP5C" ~ "Recent burn - this is the core age; but will need to update",
    TRUE ~ NA
  ))

# Here we need to make some more decisions on which tree years to use, particularly for sites that have a huge range, 
# will in general use the oldest tree unless cluster is not around that at all
# then need to check that CAFI ages were done at the tree level first too!
# I think will keep mean and max and try with both and see if it really makes a difference
rsn.diff <- rsn.trees.2 %>%
  ungroup()%>%
  group_by(Plot) %>%
  summarize(min_ty = min(mean_tree), max_ty = max(mean_tree))%>%
  mutate(diff = max_ty - min_ty) %>%
  filter(diff > 50)

ggplot(data = subset(rsn.trees.2, rsn.trees.2$Plot %in% rsn.diff$Plot), aes(mean_tree, Plot, color = SPECIES))+
  geom_point()+
  theme_bw() 


# Read in CAFI
# CAFI tree ages  - don't think need to summarize by tree; it seems that there is only one tree per sample, so should be ready to join!

cafi <- read.csv("/Users/olhajek/Desktop/RSN/RSN_proj/Data/Tree_Age_Data/Compiled/CAFI_Ages_Summarized.csv")

# Join Age summary file
glimpse(cafi)
glimpse(rsn.trees.3)

# Fix the two dataframe to merge
## CAFI
## goign to choose John Yarie and Melissa Boyd ring data when available; there's some discrepancy where it seems that Isabelle and 
## Andrew tend to have much lower ages

cafi.2 <- cafi %>%
  # choose JY and MB for 1022 adn 1088
  filter(Plot != 1022 | Windendro_Compiler != "Andrew Haverdink") %>%
  filter(Plot != 1088 | Windendro_Compiler != "Isabel Munoz") %>%
  select(-c(YearSampled, Windendro_Compiler)) %>%
  mutate(Notes = case_when(
    Plot == 1022 ~ "Removed Andrew Haverdink cores; seemed much lower",
    Plot == 1088 ~ "Removed Isabel Munoz cores; seemed much lower",
    TRUE ~ NA
  ),
  Plot = as.character(Plot)) %>%
  select(-c(mean_med, range_AGE))

rsn.3 <- rsn.trees.3 %>%
  select(-c(sd_AGE, mean_RC))

str(rsn.3)
str(cafi.2)
tree.ages <- rbind(rsn.3, cafi.2)

# Write csv - save!
write.csv(tree.ages, "/Users/olhajek/Desktop/RSN/RSN_proj/Data/Tree_Age_Data/Compiled/RSN_CAFI_Summarized_Ages.csv", row.names = FALSE)
