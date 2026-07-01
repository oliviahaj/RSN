# RSN Age checks
# Using the data from Zach to check on the ages have and create a compiled file similar to what was created with CAFI

# load libraries
library(tidyverse)

# Read in data
rsn.trees <- read.csv("/Users/olhajek/Desktop/RSN/Tree_Age_Data/Compiled/RSN_TreeAge_Compiled.csv")

# Update the dataframe for compatability with the CAFI Dataset
# need to add in AGE, rename, plot, subplot
# A few notes on tree ages
str(rsn.trees)
rsn.trees.2 <- rsn.trees %>%
  mutate(AGE = 2026 - Pith_Year) %>%
  rename(SubP = Plot, Plot = Site, YearSampled = Year_Collected, TreeYear = Pith_Year) %>%
  # Going to get rid of the Jamie for WCM1 and DCM1; for WCM2 just going to merge (113 vs. 180, different but maybe not hugely meaningfully)
  filter(Plot != "WCM1" | Source != "Jamie Found") %>%
  filter(Plot != "DCM1" | Source != "Jamie Found") %>%
  group_by(Plot) %>%
  #group_by(Plot, YearSampled, Source) %>%
  summarize(n = n(), mean_ty = mean(TreeYear, na.rm = TRUE), min_ty = min(TreeYear), max_ty = max(TreeYear),  mean_RC = mean(Ring_Count, na.rm=TRUE),
            mean_AGE = mean(AGE, na.rm = TRUE), median_AGE = median(AGE, na.rm = TRUE), min_AGE = min(AGE), max_AGE = max(AGE), sd_AGE = sd(AGE, na.rm =TRUE)) %>%
  mutate(range_AGE = max_AGE - min_AGE, mean_med = mean_AGE - median_AGE) %>%
  mutate(Notes = case_when(
    Plot == "WCM1" ~ "Removed Jamie found data; seemed to repeat data",
    Plot == "DCM1" ~ "Removed Jamie found data; huge difference from other data file; confirm",
    Plot == "FP5C" ~ "Recent burn - this is the core age; but will need to update",
    TRUE ~ NA
  ))

# Read in CAFI

cafi <- read.csv("/Users/olhajek/Desktop/RSN/Tree_Age_Data/Compiled/CAFI_Ages_Summarized.csv")

# Join Age summary file
str(cafi)
str(rsn.trees.2)

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
  Plot = as.character(Plot))

rsn.3 <- rsn.trees.2 %>%
  select(-c(sd_AGE, mean_RC))

str(rsn.3)
str(cafi.2)
tree.ages <- rbind(rsn.3, cafi.2)

# Write csv - save!
write.csv(tree.ages, "/Users/olhajek/Desktop/RSN/Tree_Age_Data/Compiled/RSN_CAFI_Summarized_Ages.csv", row.names = FALSE)
