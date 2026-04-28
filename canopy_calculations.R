## RSN DECID UNDERSTANDING

## load libraries
library(tidyverse)

# create a negate function
`%!in%` = Negate(`%in%`)

# read in data
## tree inventory
trees.3 <- read.csv('/Users/olhajek/Desktop/RSN/RSN_proj/2026 Data/QAQC_Data/Updated/RSN_Inventory_QAQC_Final_2026-04-09.csv')

## site metadata
sites <- read.csv('/Users/olhajek/Desktop/RSN/RSN_proj/2026 Data/RSN_Site/RSN_Master_Site.csv')
sites <- sites %>%
  rename(SITE = "sitename")

## filter to exclude sites not in the RSN
trees.2 <- trees.3 %>%
  filter(SITE %in% sites$SITE)

## clean up data and make sure that only trees that we want are included in this calculation
data.2 <-trees.2  %>%
  # remove salbeb
  filter(SPECIES != "SALBEB") %>%
  # remove site FP1A
  filter(SITE != "FP1A")%>%
  # set up for removing small trees
  arrange(UniqueID, Year) %>% # Ensure data is in chronological order
  group_by(UniqueID) %>%
  # keep NA DBH, DBH > 2.5, or if the DBH exceeds 2.5, it will keep subsequent observatoins
  # this code above removes 4693 observations
  filter(is.na(DBH) | DBH >= 2.49 | cumany(DBH >= 2.49)) %>%
  ungroup()


# see which trees only have an NA - this is after the seedlings have been removed, so in many cases these had an intiial measurement that was NA and corresponded to somehting like 0.7cm
na.trees <- data.2 %>%
  group_by(UniqueID) %>%
  filter(all(is.na(DBH))) #%>%
#distinct(UniqueID)

# remove the na trees only
data.25 <- data.2 %>%
  filter(UniqueID %!in% na.trees$UniqueID)


# join the data wiht the meta data
site_year <- read.csv('/Users/olhajek/Desktop/RSN/RSN_proj/2026 Data/Data_checks/site_years.csv')

data.25 <- left_join(data.25, site_year)

unique(data.25$SPECIES)

# Get rid of species that don't have allometry equatiosn from APs work
# ALSO GET RID OF DEAD/NA DBH TREES
species <- data.25 %>%
  filter(SPECIES %in% c("LARLAR", "PICEA"))

data.25.sp <- data.25 %>%
  mutate(SPECIES = ifelse(SPECIES == "PICEA", "PICMAR", SPECIES)) %>%
  filter(SPECIES != "LARLAR") %>%
  filter(!is.na(DBH))

###############################################################
#FUNCTIONS FOR TREE ALLOMETRY 
#Above DBH trees used in Alexander et al. 2012, equations from Yarie et al. 2007, function based on species and DBH (cm) measurement; output in g/tree
allom.dbh <- function(x, y) {
  treebiomass <- numeric(length(x)) 
  for (i in seq_along(x)) {
    if (x[i] == "BETNEO") {
      a <- 164.18
      b <- 2.29
    } else if (x[i] == "PICGLA") {
      a <- 96.77
      b <- 2.40
    } else if (x[i] == "PICMAR") {
      a <- 271.46
      b <- 1.84
    } else if (x[i] == "POPTRE") {
      a <- 134.10
      b <- 2.26
    } else if (x[i] == "POPBAL") {
      a <- 133.71
      b <- 2.29
    } else {
      a <- NA
      b <- NA
    }
    treebiomass[i] <- a * ((y[i])^b)
  }
  return(treebiomass) 
}


#Allometric equation for below-DBH trees - based on species and basal diameter, allometry used in Alexander et al. 2012, equations from Yarie et al. 2007; function based on species and DBH (cm) measurement; output in g/tree
allom.bd <- function(x, y) {
  treebiomass <- numeric(length(x)) 
  for (i in seq_along(x)) {
    if (x[i] == "betneo") {
      a <-  26.29
      b <- 2.68
    } else if (x[i] == "picgla") {
      a <- 53.74
      b <- 2.45
    } else if (x[i] == "picmar") {
      a <-  37.29
      b <-  2.40
    } else if (x[i] == "poptre") {
      a <-  56.83
      b <- 2.49
    } else if (x[i] == "popbal") {
      a <-  58.29
      b <- 2.44
    } else {
      a <- NA
      b <- NA
    }
    treebiomass[i] <- a * ((y[i])^b)
  }
  return(treebiomass) 
}

#Science paper allometric equation based on Jill's original seedling measurements, allometry derived from 2004 JFSP sites
allom.func.scipap <- function(x, y) {
  treebiomass <- numeric(length(x)) 
  for (i in seq_along(x)) {
    if (x[i] == "betneo") {
      a <-  2.54527
      b <- -1.2343
    } else if (x[i] == "picgla") {
      a <- 2.34579
      b <- -0.78332
    } else if (x[i] == "picmar") {
      a <-  2.61253
      b <-  -0.95462
    } else if (x[i] == "poptre") {
      a <-  2.95085
      b <- -1.57164
    } else if (x[i] == "popbal") {
      a <-  NA
      b <- NA
    } else {
      a <- NA
      b <- NA
    }
    treebiomass[i] <- 10^(((log10((y[i]*10))*a)+(b)))
  }
  return(treebiomass) 
}




#Calculate trajectories!

#Calculated by getting a mean biomass for above and below for each species, multiply by the site-level number of below and above trees, add together the biomass of below and above dbh trees of each species then divided by the area sampled to get the total biomass per msq of each vegetation type 
trajectories <- data.25.sp %>%
  #Use allometric equation to calculate biomass of individuals measured, both above and below dbh
  mutate(above_dbh_biomass = allom.dbh(SPECIES, DBH)) %>%
  #set 0 values to NA values in order to calculate the mean biomass
  mutate(above_dbh_biomass = ifelse(above_dbh_biomass == 0, NA, above_dbh_biomass)) %>%
  group_by(SITE, SPECIES, sample.rd, sample.yr, sampled_plots_period) %>%
  #for each species within each site- get the mean value for biomass for above dbh individuals and below dbh individuals. Add all below dbh individuals and all above dbh individuals for each species in each site
  summarise(mean_above_bio = mean(above_dbh_biomass, na.rm = TRUE), count_above = n()) %>%
  #multiply mean biomass * number of individuals for above and below individuals of each species, get total # of each species at the site level 
  mutate(total_above_bio = mean_above_bio * count_above, total_count = count_above) %>%
  #Get total biomass of each species at each site by getting sum of above and below biomass estimates
  mutate(total_biomass = total_above_bio) %>%
  #set spruce vs. deciduous designations by species
  mutate(veg_type = ifelse(SPECIES %in% c("BETNEO", "POPTRE", "POPBAL"), "decid", "spruce")) %>%
  group_by(SITE, sample.rd, sample.yr, sampled_plots_period, veg_type) %>%
  #get total spruce and deciduous density and biomass for each site
  summarise(site_count = sum(total_count), site_biomass = sum(total_biomass)) %>%
  #Divide the total count and biomass estimated for site by the area surveyed to get the meter squared density and biomass of each vegetation type
  #For most sites this is 14 m, for 14 quadrats measured
  #I filtered out quadrats from some sites that had clear data entry errors and have adjusted area sampled accordingly (see above for quadrats filtered out during QA/QC)
  mutate(dens_msq = site_count/sampled_plots_period) %>%
  mutate(biomass_msq = site_biomass/sampled_plots_period) %>%
  #get rid of columns you don't need 
  select(-c(site_count, site_biomass)) %>%
  #pivot data so you have a column for biomass and density for both spruce and deciduous 
  pivot_wider(names_from = veg_type, values_from = c(dens_msq, biomass_msq)) %>%
  mutate(across(starts_with("dens_msq_") | starts_with("biomass_msq_"), ~replace_na(., 0))) %>%
  #Calculate relative deciduous density, relative deciduous biomass, and deciduous index
  mutate(rel_dens_decid = dens_msq_decid/(dens_msq_decid + dens_msq_spruce), 
         rel_biomass_decid = biomass_msq_decid/(biomass_msq_decid + biomass_msq_spruce), 
         decid.index.msq = (rel_dens_decid + rel_biomass_decid)/2) %>%
  #Categorize sites based on their deciduous index value!
  mutate(canopy= ifelse((dens_msq_decid + dens_msq_spruce) < 1, "open",
                              ifelse(decid.index.msq > 0.66, "decid",
                                     ifelse(decid.index.msq > 0.33, "mixed", "spruce"))))

## Now join this with the site metadata to create a table that has the different trajectories and densities
traj.2 <- left_join(trajectories, site_year)

traj.site <- left_join(traj.2, sites, by = "SITE")


## Summary tables 
str(traj.site)
summary <- traj.site %>%
  filter(sample.yr > 2004)%>%
  group_by(sample.yr, age_class_yr, canopy) %>%
  summarize(n= n()) %>%
  pivot_wider(names_from = sample.yr,  values_from = n) 
  
summary.2 <- traj.site %>%
  filter(sample.yr > 2004)%>%
  group_by(sample.yr, age_class_yr, canopy_typ) %>%
  summarize(n= n()) %>%
  pivot_wider(names_from = sample.yr,  values_from = n) 


traj.site.map <- traj.site %>%
  group_by(SITE) %>%
  arrange(desc(sample.yr)) %>%
  slice_head(n=1)

# map
library(leaflet)
library(leaflegend)


pal <- colorFactor(palette = c("steelblue", "firebrick", "forestgreen", "goldenrod"), domain = traj.site.map$canopy)

shape_map <- c("Young" = "square", "Intermediate" = "circle", "Mature" = "triangle")

leaflet(traj.site.map) %>%
  addProviderTiles(providers$GeoportailFrance.orthos) %>% # Terrain/Satellite hybrid
  addCircleMarkers(
    lng = ~w_dd, 
    lat = ~n_dd,
    # Color by Canopy
    color = ~pal(canopy),
    fillOpacity = 0.8,
    # Change Radius based on Age Class as a visual "shape" alternative
    radius = ~ifelse(age_class_yr == "Old", 12, 6), 
    # Add a Label that appears on hover
    label = ~paste("SITE:", SITE),
    # Popup for clicking
    popup = ~paste0(
      "<b>Site: </b>", SITE, "<br>",
      "<b>Canopy: </b>", canopy, "<br>",
      "<b>Age Class: </b>", age_class_yr)) %>%
  addLegend(
    pal = pal, 
    values = ~canopy, 
    title = "Feature Type", 
    position = "bottomright"
  )
