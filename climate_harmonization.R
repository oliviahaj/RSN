# Climate compilation and cleaning
# August 27
# OLH

# load libraries
library(tidyverse)
library(leaflet)
library(viridis)

# read in the data

# monthly
month <- read.csv("/Users/olhajek/Desktop/RSN/RSN_proj/Data/climate_raw/lat_long_2012-2025MP.csv")
glimpse(month)

month.2 <- month %>%
  mutate(gs_ppt = PPT05 + PPT06 +PPT07 + PPT08, july_ppt = Tave07,  gs_temp = Tave05 + Tave06 + Tave07 + Tave08) %>%
  select(Year, ID1, gs_ppt, july_ppt, gs_temp)

# annual and seasonal
ann <- read.csv("/Users/olhajek/Desktop/RSN/RSN_proj/Data/climate_raw/lat_long_1990-2025SY.csv") 
glimpse(ann)

ann.2 <- ann %>%
  select(Year, ID1, Latitude, Longitude, Elevation, MAT, MAP, CMD, CMI, CMI_sm, CMD_sm, Tave_sp, Tave_wt, PPT_wt)

ggplot(ann.2, aes(Year, MAP, color = ID1))+
  geom_line()+
  theme_classic()+
  theme(legend.position = "none")
  
ggplot(ann.2, aes(Year, CMI_sm, color = ID1))+
  geom_line()+
  theme_classic()+
  theme(legend.position = "none")

ggplot(ann.2, aes(Year, MAT, color = ID1))+
  geom_line()+
  theme_classic()+
  theme(legend.position = "none")

# Calculate long term MAT, MAP
annual.vars <- ann %>%
  group_by(ID1, Latitude, Longitude, Elevation) %>%
  summarize(MAP = mean(MAP), MAT = mean(MAT))

# MAPs
leaflet(annual.vars) %>%
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addCircleMarkers(
    lng = ~Longitude,
    lat = ~Latitude,
    radius = 6,
    color = ~viridis(100)[cut(MAP, breaks = 100)],
    fillColor = ~viridis(100)[cut(MAP, breaks = 100)],
    fillOpacity = 0.8,
    stroke = FALSE,
    popup = ~paste0(
      "<b>Site:</b> ", ID1,
      "<br><b>MAP:</b> ", round(MAP, 1), " mm",
      "<br><b>MAT:</b> ", round(MAT, 1), " °C",
      "<br><b>Elevation:</b> ", round(Elevation, 0), " m"
    )
  )

leaflet(annual.vars) %>%
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addCircleMarkers(
    lng = ~Longitude,
    lat = ~Latitude,
    radius = 6,
    color = ~viridis(100)[cut(MAT, breaks = 100)],
    fillColor = ~viridis(100)[cut(MAT, breaks = 100)],
    fillOpacity = 0.8,
    stroke = FALSE,
    popup = ~paste0(
      "<b>Site:</b> ", ID1,
      "<br><b>MAP:</b> ", round(MAP, 1), " mm",
      "<br><b>MAT:</b> ", round(MAT, 1), " °C",
      "<br><b>Elevation:</b> ", round(Elevation, 0), " m"
    )
  )

# Join the annual and the monthly variables, then add in the long-term MAP and MAT\

clim <- left_join(month.2, ann.2)
glimpse(clim)

# try the PCA for these last few years
clim.pca <- clim %>%
  mutate(row_id = paste(ID1, Year, sep = "_")) %>%
  column_to_rownames("row_id") %>%
  select(-ID1, -Year, - Longitude, - Latitude)

pca <- prcomp(
  clim.pca,
  center = TRUE,
  scale. = TRUE
)

pca_scores <- as.data.frame(pca$x) %>%
  rownames_to_column("row_id") %>%
  separate(row_id, into = c("ID1", "Year"), sep = "_") %>%
  mutate(Year = as.integer(Year)) %>%
  left_join(clim)

ggplot(pca_scores, aes(x = PC1, y = PC2, color = as.numeric(gs_ppt))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_viridis_c() +
  theme_classic()

loadings <- as.data.frame(pca$rotation) %>%
  rownames_to_column("variable") %>%
  arrange(desc(abs(PC1)))

loadings

# individual graphs

ggplot(clim, aes(ID1, gs_ppt, color = Year))+
  geom_point()+
  theme_bw()+
  scale_color_viridis_c()
  
ggplot(clim, aes(ID1, MAP, color = Year))+
  geom_point()+
  theme_bw()+
  scale_color_viridis_c()

ggplot(clim, aes(ID1,MAT, color = Year))+
  geom_point()+
  theme_bw()+
  scale_color_viridis_c()

# Join long term data
clim.2 <- clim %>%
  rename(ann.ppt = MAP, ann.temp = MAT) %>%
  left_join(annual.vars)

# export climate data
#write.csv(clim.2, "/Users/olhajek/Desktop/RSN/RSN_proj/Data/harmonized/harmonized_climateNA.csv", row.names = F)


