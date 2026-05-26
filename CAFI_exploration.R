### CAFI Exploration - What is the coverage of the data
# How many sites
# Where are the sites
# What are the site ages?
# What are the sampling years like and how does that overlap with what we have now?

library(tidyverse)
library(leaflet)

# create a negate function
`%!in%` = Negate(`%in%`)

# read in data
cafi <- read.csv("/Users/olhajek/Desktop/RSN/CAFI_data/CAFI_Inventory_QAQC_Final.csv")
str(cafi)

# Understand if there are trees after the first sampling round that seem to big to be true
# maybe want to do a histogram to start and then see what we're working with first
cafi_new <- cafi %>%
  filter(is.na(PREVDIA))%>%
  filter(CYCLE != 1) %>%
  filter(STATUSCD != 2)%>%
  # 778 go missing
  filter(STATUSCD != 4)#%>%
  # group_by(Reconciled_Meaning)%>%
  # summarize(n = n())
# 13836 of the 16107 are ingrowth

ggplot(cafi_new, aes(DIA))+
  geom_histogram()+
  facet_wrap(~Reconciled_Meaning, scales = "free")

# everything seems to make sense but not all of the NAs do so will look at that more
# there are 1570 observations - any patterns; it seems like a lot of these would be ingrowth
# a lot of these have measured correctly in there..
cafi_na <- cafi_new %>%
  filter(is.na(Reconciled_Meaning))

# should redo these analyses and look at this once I have a sense of which sites should be included as sites in interior :) 

# What is the forest type for these different sites and what are their ages? Where is this info?
# This will be good info to show on the visualizations later
locs <- read.csv("/Users/olhajek/Desktop/RSN/CAFI_data/452_CAFI_LOCATION_FUZZED_v20241220.csv")
str(locs)
locs.avg <- locs %>%
  #a weird value that messes up the whole 1098 plots
  filter(OID_ !=290)%>%
  # plot1180 has lat and long backwards on one subplot
  mutate(LAT = ifelse(OID_ == 534, LONG, LAT),
         LONG = ifelse(OID_ == 534, -149.2414022, LONG))%>%
  group_by(PLOT) %>%
  summarize(lat = mean(LAT), long= mean(LONG), n = n()) %>%
  # get rid of a plots that don't have teh right lat/long - not sure what's going on with these
  # check these plots in the data - were they discontinued?
  filter(PLOT %!in% c(1206, 1040, 1026))
# note that there are two sites that only have one subplot in this

leaflet(locs.avg) %>%
  addProviderTiles("Esri.WorldTopoMap") %>%
  addCircleMarkers(
    lng = ~long,
    lat = ~lat,
    popup = ~as.character(PLOT)
  )


locs.zach <- read.csv("/Users/olhajek/Desktop/RSN/CAFI_data/cafi_loc_deciduous.csv")
leaflet(locs.zach) %>%
  addProviderTiles("Esri.WorldTopoMap") %>%
  addCircleMarkers(
    lng = ~Longitude,
    lat = ~Latitude,
    popup = ~as.character(site)
  )

# what is missing in the overlap
overlap <- locs.avg %>%
  filter(PLOT %!in% locs.zach$site)

check <- locs.avg %>%
  filter(PLOT %in% c(1141, 1149))

leaflet(check) %>%
  addProviderTiles("Esri.WorldTopoMap") %>%
  addCircleMarkers(
    lng = ~long,
    lat = ~lat,
    popup = ~as.character(PLOT)
  )


# which sites to include!!
locs.interior <- locs.avg %>%
  filter(lat > 62.69353) %>%
  filter(!(lat < 63.15401 & long < -143.7223))

leaflet(locs.interior) %>%
  addProviderTiles("Esri.WorldTopoMap") %>%
  addCircleMarkers(
    lng = ~long,
    lat = ~lat,
    popup = ~as.character(PLOT)
  )

# looks like there are 78 plots in interior alaska 
# now filter the CAFI dataset to just this 
cafi.interior <- cafi %>%
  filter(PLOT %in% locs.interior$PLOT)

# look at the number of "new" trees
cafi_new <- cafi.interior %>%
  filter(is.na(PREVDIA))%>%
  filter(CYCLE != 1) %>%
  filter(STATUSCD != 2)%>%
  # 778 go missing
  filter(STATUSCD != 4) #%>%
  # group_by(Reconciled_Meaning)%>%
  # summarize(n = n())
# 6105 of the 6872 are ingrowth
# 1570? NAs

ggplot(cafi_new, aes(DIA))+
  geom_histogram()+
  facet_wrap(~Reconciled_Meaning, scales = "free")

cafi.na <- cafi_new %>%
  filter(is.na(Reconciled_Meaning))

# Look at each site over time with the trees classed by alive, dead, etc. 
unique(cafi.interior$STATUSCD)
data_summary <- cafi.interior %>%
  # make a column for whether or not tree is above or below 2.5 cm
  mutate(over2.5 = ifelse(DIA>2.5, 1, 0)) %>%
  # make a column that indicates, dead, live, or seedling less than 2.5 cm
  mutate(tree_status = case_when(
    over2.5 == 1 & STATUSCD== 1 ~ as.character("Live"), 
    over2.5 == 0 & STATUSCD ==1 ~ as.character("Seedling"),
    STATUSCD == 2 ~ as.character("Dead"),
    STATUSCD == 4 ~ as.character("Missing"), 
    STATUSCD == 3 ~ as.character("Harvested"),
    TRUE ~ "TBD"
  )) %>%
  group_by(PLOT, INVYR, tree_status) %>%
  summarize(no_trees = n_distinct(TREE_ID)) %>%
  pivot_wider(names_from = tree_status, values_from = no_trees, values_fill = 0) %>%
  ungroup()%>%
  group_by(PLOT) %>%
  arrange(INVYR, .by_group = TRUE) %>%
  mutate(live_diff = Live - lag(Live))

# 2. Get a list of all sites
plot_ready <- data_summary %>%
  # Select only the columns needed for the graph
  select(PLOT, INVYR, Live, Dead, Seedling, Missing, TBD) %>%
  # Pivot so 'Live' and 'Dead' are in one column
  pivot_longer(cols = c(Live, Dead, Seedling, Missing, TBD), 
               names_to = "Status", 
               values_to = "Count") %>%
  # Clean up labels: "no_trees_Live" becomes "Live"
  mutate(Status = str_remove(Status, "no_trees_"))
site_list <- unique(plot_ready$PLOT)

for (s in site_list) {
  p <- plot_ready %>%
    filter(PLOT == s) %>%
    ggplot(aes(x = INVYR, y = Count, color = Status, group = Status)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = c("Live" = "#2ecc71", "Dead" = "#e74c3c", "Seedling" = "skyblue",
                                  "Missing" = "grey", "TBD" = "purple")) +
    labs(title = paste("Tree Density Trends - Site:", s),
         subtitle = "Comparing Live vs. Dead stems over time",
         x = "Year", y = "Number of Trees") +
    theme_minimal()
  
  print(p) # This "prints" the plot as a new page in the PDF
}


# will want to do this again once i actually get the plots and stuff with multiple measurements, etc

# What do the sampling years look like and how does that overlap with what we have now?
# what info do we have on age class and canopy class
locs.zach <- locs.zach %>%
  rename(PLOT = "site")
cafi.interior <- left_join(cafi.interior, locs.zach)

tbd <- cafi.interior %>%
  filter(is.na(Class)) %>%
  select(PLOT)%>%
  distinct()

str(cafi.interior)
cafi.summary <- cafi.interior %>%
  mutate(Class = ifelse(is.na(Class), "UNK", Class)) %>%
  group_by(PLOT, Class) %>%
  summarize(first_year = min(INVYR), last_year = max(INVYR), num_inv = length(unique(CYCLE)), age = mean(StdAge)) %>%
  mutate(age2 = ifelse(PLOT %in% c(1007, 1008, 1011, 1013, 1014, 1018, 1019, 1020, 1021, 1022, 1080, 1081,
                                   1082, 1083, 1084, 1090, 1091, 1101, 1102, 1103, 1104, 1146, 1187), "Coming", "Not"))

cafi.samples <- cafi.interior %>%
  mutate(Class = ifelse(is.na(Class), "UNK", Class))  %>%
  select(PLOT, Class, INVYR)
unique(cafi.samples$Class)

ggplot(data = subset(cafi.samples, cafi.samples$Class == "Decid"), aes(as.factor(PLOT), INVYR, color = Class)) + 
  geom_point() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8))

ggplot(data = subset(cafi.samples, cafi.samples$Class == "Spruce"), aes(as.factor(PLOT), INVYR, color = Class)) + 
  geom_point() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8))

ggplot(data = subset(cafi.samples, cafi.samples$Class == "Mixed"), aes(as.factor(PLOT), INVYR, color = Class)) + 
  geom_point() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8))

ggplot(data = subset(cafi.samples, cafi.samples$Class == "UNK"), aes(as.factor(PLOT), INVYR, color = Class)) + 
  geom_point() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8))

# Which sites have a measurement in 2020
cafi.2020s <- cafi.interior %>%
  filter(INVYR>2019) %>%
  select(PLOT)%>%
  distinct()

cafi.subset <- cafi.interior%>%
  filter(PLOT %in% cafi.2020s$PLOT)

ggplot(data = cafi.subset, aes(as.factor(PLOT), INVYR, color = Class)) + 
  geom_point() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8))

ggplot(data = cafi.subset, aes(as.factor(PLOT), INVYR, color = Class)) + 
  geom_point() +
  theme_bw()+
  ylim(2011, 2025)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8))

