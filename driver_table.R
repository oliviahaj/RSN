# Driver file creation
# OLH
# Oct 27

# load libraries 
library(tidyverse)

# read in data
## climate
climate <- read.csv("/Users/olhajek/Desktop/RSN/RSN_proj/Data/harmonized/harmonized_climateNA.csv")
sampling <- read.csv("/Users/olhajek/Desktop/RSN/RSN_proj/Data/harmonized/continuous_tree_ages.csv")
site <- read.csv('/Users/olhajek/Desktop/RSN/RSN_proj/Data/2026 Data/RSN_Site/RSN_Master_Site.csv')
glimpse(climate)

# sampling summary 
glimpse(sampling)
sampling2 <- sampling %>%
  select(c(PLOT, Year, sample.pd)) 

# Use the smapling summary with the climate to calculate the mean climate over the different windows
periods <- sampling2 %>%
  arrange(PLOT, Year) %>%
  group_by(PLOT) %>%
  mutate(end   = Year,
         start = coalesce(lag(Year) + 1L, 2000)) |>
  ungroup() %>%
  select(PLOT, sample.pd, start, end) %>%
  filter(sample.pd != 1) %>%
  mutate(start = case_when(
    sample.pd == 3 & start == 2000 ~ 2018,
    sample.pd == 2 & start == 2000 & PLOT == "FP5A" ~ 2010,
    sample.pd == 2 & start == 2000 & PLOT == "FP5D" ~ 2010,
    sample.pd == 2 & start == 2000 ~ 2015,
    TRUE ~ start
  )) %>%
  mutate(end2 = end -1) %>%
  select(-end)

glimpse(climate)
clim_means <- periods |>
  inner_join(climate, join_by(PLOT == ID1, start <= Year, end2 >= Year)) |>
  group_by(PLOT, sample.pd, start, end2) |>
  summarise(n_yrs = n(), max_julytemp = max(july_maxtemp), max_gs_temp = max(gs_temp), 
            max_CMI_sm = max(CMI_sm),
            across(gs_ppt:MAT & !c(Latitude, Longitude, Elevation), mean),
            .groups = "drop")


# Join with the site data
glimpse(clim_means)
glimpse(site)

unique(site$moist)
#clean up site
site2 <- site %>%
  select(sitename, slope, aspect, elevation, topo, moist, depth.mean, depth.std, 
         depth.rep, ph, n_dd, w_dd) %>%
  rename(PLOT = sitename) %>%
  mutate(moist_class = case_when(
    moist == "xeric" ~ 1, 
    moist == "subxeric" ~ 2,
    moist == "subxeric-mesic" ~ 3,
    moist == "mesic" ~ 4,
    moist == "mesic-subhygric" ~ 5,
    TRUE ~ 6
  ))

site_data <- left_join(clim_means, site2, by = "PLOT")

## SITE CONDITIONS PCA w/o ph

# Make some figures and PCA of site conditions...
glimpse(site_data)

site_pca <- site_data |>
  distinct(PLOT, slope, aspect, elevation, moist_class) %>%
  mutate(northness = if_else(slope == 0, 0, cos(aspect * pi / 180)),
         eastness  = if_else(slope == 0, 0, sin(aspect * pi / 180)))

# try the PCA for these last few years
site_pca2 <- site_pca %>%
  column_to_rownames("PLOT") %>%
  select(-aspect)

pca <- prcomp(
  site_pca2,
  center = TRUE,
  scale. = TRUE
)

pca_scores <- as.data.frame(pca$x) %>%
  rownames_to_column("PLOT") %>%
  left_join(site_pca)

glimpse(pca_scores)
ggplot(pca_scores, aes(x = PC1, y = PC2, color = as.numeric(elevation))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_viridis_c() +
  theme_classic()

ggplot(pca_scores, aes(x = PC1, y = PC2, color = moist_class)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_viridis_c() +
  theme_classic()

ggplot(pca_scores, aes(x = PC1, y = PC2, color = slope)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_viridis_c() +
  theme_classic()

ggplot(pca_scores, aes(x = PC1, y = PC2, color = northness)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_viridis_c() +
  theme_classic()

ggplot(pca_scores, aes(x = PC1, y = PC2, color = eastness)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_viridis_c() +
  theme_classic()

loadings <- as.data.frame(pca$rotation) %>%
  rownames_to_column("variable") %>%
  arrange(desc(abs(PC1)))

loadings
summary(pca)

### NOT GOING TO INCLUDE PH
## PCA w/ PH and subset of sites
# site_pca <- site_data |>
#   distinct(PLOT, slope, aspect, elevation, moist_class, ph) %>%
#   mutate(northness = if_else(slope == 0, 0, cos(aspect * pi / 180)),
#          eastness  = if_else(slope == 0, 0, sin(aspect * pi / 180))) %>%
#   filter(!is.na(ph))
# 
# # try the PCA for these last few years
# site_pca2 <- site_pca %>%
#   column_to_rownames("PLOT") %>%
#   select(-c(aspect))
# 
# pca <- prcomp(
#   site_pca2,
#   center = TRUE,
#   scale. = TRUE
# )
# 
# pca_scores <- as.data.frame(pca$x) %>%
#   rownames_to_column("PLOT") %>%
#   left_join(site_pca)
# 
# glimpse(pca_scores)
# ggplot(pca_scores, aes(x = PC1, y = PC2, color = as.numeric(elevation))) +
#   geom_point(size = 3, alpha = 0.7) +
#   scale_color_viridis_c() +
#   theme_classic()
# 
# loadings <- as.data.frame(pca$rotation) %>%
#   rownames_to_column("variable") %>%
#   arrange(desc(abs(PC1)))
# 
# loadings
# summary(pca)

# Join the PCA 
glimpse(pca_scores)

pca_join <- pca_scores %>%
  select(c(PLOT, PC1, PC2)) %>%
  rename(topo_moist_PC1 = PC1, 
         cold_PC2 = PC2)

site_data2 <- site_data %>%
  left_join(pca_join)

## Explore the climate periods
glimpse(site_data2)

ggplot(site_data2, aes(max_CMI_sm, max_gs_temp, group = PLOT, color = n_dd))+
  geom_point(, size = 2)+
  geom_line()+
  theme_bw() +
  facet_wrap(~sample.pd)


glimpse(climate)
ggplot(climate, aes(Year, july_maxtemp , color = ID1))+
  geom_line()+
  theme_bw()

ggplot(climate, aes(Year, july_maxtemp , color = ID1))+
  geom_line()+
  theme_bw()

ggplot(climate, aes(Year, CMI_sm, color = ID1))+
  geom_line()+
  theme_bw()

ggplot(climate, aes(Year, july_ppt , color = ID1))+
  geom_line()+
  theme_bw()

ggplot(climate, aes(Year, ann.temp , color = ID1))+
  geom_line()+
  theme_bw()

ggplot(climate, aes(Year, ann.ppt , color = ID1))+
  geom_line()+
  theme_bw()

ggplot(climate, aes(Year, CMI , color = ID1))+
  geom_line()+
  theme_bw()

#write.csv(site_data2, '/Users/olhajek/Desktop/RSN/RSN_proj/Data/2026 Data/RSN_Site/RSN_Climate_Site_characteristcs.csv', row.names = F)
