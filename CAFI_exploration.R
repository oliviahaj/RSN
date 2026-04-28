### CAFI Exploration - What is the coverage of the data

library(tidyverse)

# read in data
cafi <- read.csv("/Users/olhajek/Desktop/RSN/CAFI_data/CAFI_Inventory_Meta-Added.csv")
cafi25 <- read.csv("/Users/olhajek/Desktop/RSN/CAFI_data/cafi_2025_sites.csv")
cafi_meta <- read.csv("/Users/olhajek/Desktop/RSN/CAFI_data/cafi_loc_deciduous.csv")

# Extract site-year 
cafi_sample_yrs <- cafi %>%
  select(PLOT, INVYR) %>%
  distinct()

caf25 <- cafi25 %>%
  rename(PLOT2 = "Site") %>%
  mutate(PLOT = substr(PLOT2, nchar(PLOT2) - 3, nchar(PLOT2)), 
         INVYR = 2025) %>%
  select(-PLOT2)

cafi_samples <- rbind(cafi_sample_yrs, caf25)

cafi_meta <- cafi_meta %>%
  rename(PLOT = "site")%>%
  mutate(PLOT = as.character(PLOT)) %>%
  mutate(Class = factor(Class, 
                            levels = c("Decid", "Mixed", "Spruce", "NA")))

cafi_samples.meta <- left_join(cafi_samples, cafi_meta)

## Visualize - Overall - OMG!
ggplot(cafi_samples.meta, aes(PLOT,INVYR, color = Class))+
  geom_point()+
  theme_bw()


cafi_sample.meta2 <- cafi_samples.meta 

ggplot(data = subset(cafi_sample.meta2, cafi_sample.meta2$Class == "Decid"), aes(PLOT,INVYR, color = Class))+
  geom_point(size = 3)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14)) +
  ggtitle("Decid Stands")

ggplot(data = subset(cafi_sample.meta2, cafi_sample.meta2$Class == "Mixed"), aes(PLOT,INVYR, color = Class))+
  geom_point(size = 3)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14)) +
  ggtitle("Mixed Stands")

ggplot(data = subset(cafi_sample.meta2, cafi_sample.meta2$Class == "Spruce"), aes(PLOT,INVYR, color = Class))+
  geom_point(size = 3)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14)) +
  ggtitle("Spruce Stands")

ggplot(data = subset(cafi_sample.meta2, is.na(cafi_sample.meta2$Class)), aes(PLOT,INVYR, color = Class))+
  geom_point(size = 3)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14)) +
  ggtitle("Unclassified Stands")

# Table to summarize
summary <- cafi_samples.meta %>%
  filter(INVYR > 2010) %>%
  group_by(INVYR, Class) %>%
  summarize(n= n()) %>%
  pivot_wider(names_from = INVYR,  values_from = n) 


