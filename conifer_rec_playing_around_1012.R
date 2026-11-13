# Doing the biomass accumulation model
# Nov 9, 2026

# load libraries
library(tidyverse)
library(glmmTMB)
library(splines)
library(DHARMa)

# create a negate function
`%!in%` = Negate(`%in%`)

# Read in the data
## Inventory data
trees <- read.csv("/Users/olhajek/Desktop/RSN/RSN_proj/Data/harmonized/RSN_CAFI_merged.csv")

## Tree ages## Tree agesdb_has_table()
ages <- read.csv("/Users/olhajek/Desktop/RSN/RSN_proj/Data/harmonized/continuous_tree_ages.csv") 

# some data cleaning
qaqc <- trees %>%
  filter(QAQC_Action %in% c("Flag_LastObsExtremeDip",
                            # going to keep this code in right now, but not quite how this code is fully being populated
                            # "Flag_NewGrowth_Maybe?", 
                            "Flag_LastObs_SpeciesGrowthSpike"))

qaqctrees <- trees %>%
  filter(UniqueID %in% qaqc$UniqueID)

trees.2 <- trees %>%
  filter(UniqueID %!in% qaqc$UniqueID) %>%
  # this is a single observation that is in the wrong year and doesn't seem quite right
  filter(UniqueID != "WCM2_14_2939")

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
      # Making PICEA use the PICGLA allometry equations b/c PICEA is only at one site and no PICMAR recorded there, only PICMAR
    } else if (x[i] == "PICEA") {
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
    } else if (x[i] == "LARLAR") {
      a <- 10^2.648
      b <- 0.715
    } else {
      a <- NA
      b <- NA
    }
    treebiomass[i] <- a * ((y[i])^b)
  }
  return(treebiomass) 
}


# decid vs. conifer set up
decid <- c("POPTRE", "BETNEO", "POPBAL")

# create a dataframe that has the deciduous index, plot and age
# for sites that only have one interval, need
# live trees
tree.summary.species <- trees.2 %>%
  mutate(above_dbh_biomass = allom.dbh(SPECIES, DBH)) %>%
  left_join(ages) %>%
  mutate(func_gp = ifelse(SPECIES %in% decid, "Decid_Broad", "Conifer")) %>%
  filter(DA!=0) %>%
  mutate(area.sampled = Plot.Area.M2 * sampled_plots_period, 
         BA = 0.00007854 * (DBH)^2) %>%
  group_by(func_gp,Year, PLOT, area.sampled) %>%
  summarize(n = n(), mean_dbh = mean(DBH, na.rm = T), se_dbh = sd(DBH, na.rm = T)/sqrt(n()), 
            ba = sum(BA, na.rm = T), mean_above_bio = sum(above_dbh_biomass, na.rm = TRUE)) %>%  
  ungroup() %>%
  mutate(across(where(is.numeric), \(x) replace_na(x, 0))) %>%
  mutate(dens = n/area.sampled, ba_m2 = ba/area.sampled, agb_m2 = mean_above_bio/area.sampled) %>%
  pivot_wider(id_cols = c(Year, PLOT, area.sampled),
              names_from = func_gp,
              values_from = c(dens, mean_dbh, se_dbh, ba_m2, agb_m2),
              names_glue = "{func_gp}_{.value}") %>%
  mutate(across(ends_with(c("_dens", "_ba_m2", "mean_dbh", "se_dbh", "agb_m2")), ~replace_na(.x, 0))) 

tree.summary <- trees.2 %>%
  mutate(above_dbh_biomass = allom.dbh(SPECIES, DBH)) %>%
  left_join(ages) %>%
  mutate(func_gp = ifelse(SPECIES %in% decid, "Decid_Broad", "Conifer")) %>%
  filter(DA!=0) %>%
  mutate(area.sampled = Plot.Area.M2 * sampled_plots_period, 
         BA = 0.00007854 * (DBH)^2) %>%
  group_by(Year, PLOT, area.sampled) %>%
  summarize(n = n(), mean_dbh = mean(DBH, na.rm = T), se_dbh = sd(DBH, na.rm = T)/sqrt(n()), 
            ba = sum(BA, na.rm = T), mean_above_bio = sum(above_dbh_biomass, na.rm = TRUE)) %>%  
  ungroup() %>%
  mutate(across(where(is.numeric), \(x) replace_na(x, 0))) %>%
  mutate(dens = n/area.sampled, ba_m2 = ba/area.sampled, agb_m2 = mean_above_bio/area.sampled) %>%
  select(-c(n, ba, mean_above_bio))

tree.sum <- left_join(tree.summary, tree.summary.species) %>%
  mutate(across(where(is.numeric), \(x) replace_na(x, 0))) %>%
  #Calculate relative deciduous density, relative deciduous biomass, and deciduous index
  mutate(rel_dens_decid = Decid_Broad_dens/dens, 
         rel_biomass_decid = Decid_Broad_agb_m2/agb_m2, 
         decid.index.msq = (rel_dens_decid + rel_biomass_decid)/2) %>%
  left_join(ages) %>%
  mutate(age2 = Year - age)

# now I need to trim down the dataset and then I need to add in the 0 sampling efforts
glimpse(tree.sum)
tree.sum2 <- tree.sum %>%
  select(c(Network, Year, PLOT, decid.index.msq, sample.pd, age2))

# look at decid index and age
ggplot(tree.sum2, aes(age2, decid.index.msq, color = Network))+
  geom_point()+
  theme_bw() 

#### Need to update so that all of the intervals have a pre....
# some young sites that don't have a first interval with 0 trees....
tree.ct <- tree.sum2 %>%
  group_by(PLOT) %>%
  summarize(n = n(), min_sp = min(sample.pd)) %>%
  filter(n<3) %>%
  filter(PLOT %!in% c("BDI5", "BDI6", "FP5C", "FP5A", "FP5D")) %>%
  filter(!str_detect(PLOT, "M"))

tree.miss <- tree.sum2 %>%
  filter(PLOT %in% tree.ct$PLOT) %>%
  mutate(plot.sr2 = paste(PLOT, sample.pd, sep = "_"))

tree.miss.join <- tree.miss %>%
  group_by(PLOT) %>%
  slice_min(Year, n = 1) %>%
  mutate(Year2 = Year) %>%
  select(-c(Year,sample.pd, plot.sr2, Network))

plots <- tree.miss %>%
  select(PLOT)
sample.rd <- tree.miss %>%
  select(sample.pd) %>%
  unique()

plot.sr <- cross_join(plots, sample.rd) %>%
  mutate(plot.sr2 = paste(PLOT, sample.pd, sep = "_")) %>%
  filter(plot.sr2 %!in% tree.miss$plot.sr2) %>%
  mutate(Year = case_when(
    sample.pd == 1 ~ 2014, 
    sample.pd == 3 ~ 2024, 
    TRUE ~ 2018
  )) %>%
  select(-plot.sr2) %>%
  left_join(tree.miss.join) %>%
  mutate(year_diff = Year2-Year, age3 = age2 - year_diff) %>%
  select(-c(year_diff, Year2, age2)) %>%
  rename(age2 = age3) %>%
  mutate(Network = "RSN")

# join the two
tree.meta <- rbind(tree.sum2, plot.sr)

drivers <- read.csv('/Users/olhajek/Desktop/RSN/RSN_proj/Data/2026 Data/RSN_Site/RSN_Climate_Site_characteristcs.csv') 

pc1 <- drivers %>%
  select(c(PLOT, topo_moist_PC1, cold_PC2)) %>%
  unique()

tree.meta2 <- left_join(tree.meta, pc1) %>%
  filter(PLOT %!in% c("BDI5", "BDI6", "FP5C", "FP5A", "FP5D")) %>%
  distinct()


#########################
# tree and biomass
glimpse(trees.2)
tree.subplot <- trees.2 %>%
  mutate(above_dbh_biomass = allom.dbh(SPECIES, DBH),
         BA = 0.00007854 * (DBH)^2) %>%
  filter(PLOT %!in% c("BDI5", "BDI6", "FP5C", "FP5A", "FP5D")) %>%
  filter(Network == "RSN") %>%
  filter(DA == 1) %>%
  group_by(Year, PLOT, SUBP) %>%
  summarize(biomass = sum(above_dbh_biomass, na.rm = T), ba = sum(BA, na.rm = T), n=n()) %>%
  ungroup()%>%
  mutate(biomass_m2_sp = biomass/100, ba_m2_sp = ba/100, dens_sp = n/100) %>%
  select(-c(biomass, ba, n))

subplot_traj <- trees.2 %>%
  mutate(above_dbh_biomass = allom.dbh(SPECIES, DBH),
         BA = 0.00007854 * (DBH)^2) %>%
  mutate(func_gp = ifelse(SPECIES %in% decid, "Decid_Broad", "Conifer")) %>%
  filter(PLOT %!in% c("BDI5", "BDI6", "FP5C", "FP5A", "FP5D")) %>%
  filter(Network == "RSN") %>%
  filter(DA == 1) %>%
  group_by(Year, PLOT, SUBP, func_gp) %>%
  summarize(biomass = sum(above_dbh_biomass), ba = sum(BA), n=n()) %>%
  ungroup()%>%
  mutate(biomass_m2_sp = biomass/100, ba_m2_sp = ba/100, dens_sp = n/100) %>%
  select(-c(biomass, ba, n)) %>%
  pivot_wider(id_cols = c(Year, PLOT, SUBP),
              names_from = func_gp,
              values_from = c(dens_sp, ba_m2_sp, biomass_m2_sp),
              names_glue = "{func_gp}_{.value}") %>%
  mutate(across(ends_with(c("_dens_sp", "_ba_m2_sp", "biomass_m2_sp")), ~replace_na(.x, 0))) %>% 
  mutate(rel_dens_decid = Decid_Broad_dens_sp/(Decid_Broad_dens_sp + Conifer_dens_sp), 
         rel_biomass_decid = Decid_Broad_biomass_m2_sp/(Decid_Broad_biomass_m2_sp + Conifer_biomass_m2_sp), 
         decid.index.msq.sp = (rel_dens_decid + rel_biomass_decid)/2) %>%
  select(-c(rel_dens_decid, rel_biomass_decid)) %>%
  mutate(decid.index.msq.sp = replace_na(decid.index.msq.sp, 0))%>%
  select(c(Year, PLOT, SUBP, decid.index.msq.sp))

tree.subplot2 <- left_join(tree.subplot, subplot_traj)

ggplot(tree.subplot, aes(PLOT, dens))+
  geom_point()+
  theme_bw()+
  facet_wrap(~Year)
ggplot(tree.subplot, aes(PLOT,biomass_m2))+
  geom_point()+
  theme_bw()+
  facet_wrap(~Year)

# join subplot
tree.combo <- left_join(tree.subplot2, tree.sum)

# Visualize
glimpse(tree.combo)

ggplot(tree.combo, aes(age2, dens))+
  geom_point(aes(age2, dens_sp, color = decid.index.msq), alpha = 0.2)+
  geom_line(aes(group = PLOT)) +
  geom_point(aes(color = decid.index.msq), size =2 ) +
  theme_bw() +
  xlim(0, 250)+
  ylab("Density m2")+
  xlab("Time since fire")

ggplot(tree.combo, aes(age2, agb_m2))+
  geom_point(aes(age2, biomass_m2_sp, color = decid.index.msq), alpha = 0.2)+
  geom_line(aes(group = PLOT)) +
  geom_point(aes(color = decid.index.msq), size =2 ) +
  theme_bw() +
  xlim(0, 250)+
  ylab("AGB m2")+
  xlab("Time since fire")

library(ggnewscale)

decid_scale <- function(name) {
  scale_color_gradientn(
    colours = c("#1b4d3e", "#7fa37a", "#e8d9a0", "#d4a017", "#b5651d"),
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c("Conifer", "Mixed", "Deciduous"),
    name = name
  )
}

ggplot(tree.combo, aes(age2, agb_m2)) +
  geom_point(aes(y = biomass_m2_sp, color = decid.index.msq.sp), alpha = 0.2) +
  decid_scale("Species") +
  new_scale_color() +
  geom_line(aes(group = PLOT), color = "grey60") +
  geom_point(aes(color = decid.index.msq), size = 2) +
  decid_scale("Plot") +
  coord_cartesian(xlim = c(0, 250)) +
  labs(x = "Time since fire", y = "AGB m2") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

pal <- c("#1b4d3e", "#7fa37a", "#e8d9a0", "#d4a017", "#b5651d")

ggplot(tree.combo, aes(age2, agb_m2)) +
  geom_point(aes(y = biomass_m2_sp, color = decid.index.msq.sp),
             shape = 1, size = 1, alpha = 0.5) +
  geom_line(aes(group = PLOT), color = "grey40", linewidth = 0.4) +
  geom_point(aes(fill = decid.index.msq),
             shape = 21, size = 2.8, color = "white", stroke = 0.4) +
  scale_color_gradientn(colours = pal, limits = c(0, 1),
                        breaks = c(0, 0.5, 1),
                        labels = c("Conifer", "Mixed", "Deciduous"),
                        name = "Species") +
  scale_fill_gradientn(colours = pal, limits = c(0, 1),
                       breaks = c(0, 0.5, 1),
                       labels = c("Conifer", "Mixed", "Deciduous"),
                       name = "Plot") +
  coord_cartesian(xlim = c(0, 250)) +
  labs(x = "Time since fire", y = "AGB m2") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())


ggplot(tree.combo, aes(age2, agb_m2)) +
  # geom_point(aes(y = biomass_m2_sp, color = decid.index.msq.sp),
  #            shape = 1, size = 1, alpha = 0.5) +
  geom_line(aes(group = PLOT), color = "grey40", linewidth = 0.4) +
  geom_point(aes(fill = decid.index.msq),
             shape = 21, size = 2.8, color = "white", stroke = 0.4) +
  # scale_color_gradientn(colours = pal, limits = c(0, 1),
  #                       breaks = c(0, 0.5, 1),
  #                       labels = c("Conifer", "Mixed", "Deciduous"),
  #                       name = "Species") +
  scale_fill_gradientn(colours = pal, limits = c(0, 1),
                       breaks = c(0, 0.5, 1),
                       labels = c("Conifer", "Mixed", "Deciduous"),
                       name = "Plot") +
  coord_cartesian(xlim = c(10, 65)) +
  labs(x = "Time since fire", y = "AGB m2") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())



# check site quality and AGB
tree.combo <-left_join(tree.combo, pc1)

pal_moist <- c("#08306b", "#2c7fb8", "#7fcdbb", "#edf8b1", "#d8b365", "#8c510a")


ggplot(tree.combo, aes(age2, agb_m2)) +
  # geom_point(aes(y = biomass_m2_sp, color = decid.index.msq.sp),
  #            shape = 1, size = 1, alpha = 0.5) +
  geom_line(aes(group = PLOT), color = "grey40", linewidth = 0.4) +
  geom_point(aes(fill = topo_moist_PC1),
             shape = 21, size = 2.8, color = "white", stroke = 0.4) +
  # scale_color_gradientn(colours = pal, limits = c(0, 1),
  #                       breaks = c(0, 0.5, 1),
  #                       labels = c("Conifer", "Mixed", "Deciduous"),
  #                       name = "Species") +
  scale_fill_gradientn(colours = pal_moist, 
                       breaks = c(-1.5, 0.5, 2.5),
                       labels = c("Wetter/Flatter", "Mesic", "Drier/Steeper"),
                       name = "Plot") +
  coord_cartesian(xlim = c(10, 150)) +
  labs(x = "Time since fire", y = "AGB m2") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

ggplot(tree.combo, aes(age2, agb_m2)) +
  # geom_point(aes(y = biomass_m2_sp, color = decid.index.msq.sp),
  #            shape = 1, size = 1, alpha = 0.5) +
  geom_line(aes(group = PLOT), color = "grey40", linewidth = 0.4) +
  geom_point(aes(fill = cold_PC2),
             shape = 21, size = 2.8, color = "white", stroke = 0.4) +
  # scale_color_gradientn(colours = pal, limits = c(0, 1),
  #                       breaks = c(0, 0.5, 1),
  #                       labels = c("Conifer", "Mixed", "Deciduous"),
  #                       name = "Species") +
  scale_fill_gradientn(colours = pal_moist, 
                       breaks = c(-1.5, 0.5, 2.5),
                       labels = c("Wetter/Flatter", "Mesic", "Drier/Steeper"),
                       name = "Plot") +
  coord_cartesian(xlim = c(10, 150)) +
  labs(x = "Time since fire", y = "AGB m2") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

## DENSITY
ggplot(tree.combo, aes(age2, dens)) +
  geom_point(aes(y = dens_sp, color = decid.index.msq.sp),
             shape = 1, size = 1, alpha = 0.5) +
  geom_line(aes(group = PLOT), color = "grey40", linewidth = 0.4) +
  geom_point(aes(fill = decid.index.msq),
             shape = 21, size = 2.8, color = "white", stroke = 0.4) +
  scale_color_gradientn(colours = pal, limits = c(0, 1),
                        breaks = c(0, 0.5, 1),
                        labels = c("Conifer", "Mixed", "Deciduous"),
                        name = "Species") +
  scale_fill_gradientn(colours = pal, limits = c(0, 1),
                       breaks = c(0, 0.5, 1),
                       labels = c("Conifer", "Mixed", "Deciduous"),
                       name = "Plot") +
  coord_cartesian(xlim = c(0, 250)) +
  labs(x = "Time since fire", y = "Density m2") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())


ggplot(tree.combo, aes(age2, dens)) +
  # geom_point(aes(y = biomass_m2_sp, color = decid.index.msq.sp),
  #            shape = 1, size = 1, alpha = 0.5) +
  geom_line(aes(group = PLOT), color = "grey40", linewidth = 0.4) +
  geom_point(aes(fill = decid.index.msq),
             shape = 21, size = 2.8, color = "white", stroke = 0.4) +
  # scale_color_gradientn(colours = pal, limits = c(0, 1),
  #                       breaks = c(0, 0.5, 1),
  #                       labels = c("Conifer", "Mixed", "Deciduous"),
  #                       name = "Species") +
  scale_fill_gradientn(colours = pal, limits = c(0, 1),
                       breaks = c(0, 0.5, 1),
                       labels = c("Conifer", "Mixed", "Deciduous"),
                       name = "Plot") +
  coord_cartesian(xlim = c(10, 65)) +
  labs(x = "Time since fire", y = "Density m2") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

## Dying and recruitment

# prep dataframe
tree.demo.subp <- trees.2 %>%
  mutate(func_gp = ifelse(SPECIES %in% decid, "Decid_Broad", "Conifer")) %>%
  mutate(above_dbh_biomass = allom.dbh(SPECIES, DBH),
         BA = 0.00007854 * (DBH)^2) %>%
  filter(PLOT %!in% c("BDI5", "BDI6", "FP5C", "FP5A", "FP5D")) %>%
  filter(Network == "RSN") %>%
  mutate(status = ifelse(DA == 1, "live", "dead"))

subplot_traj 


# ---- 1. Census intervals per plot ---------------------------------------------
intervals <- tree.demo.subp |>
  select(PLOT, SUBP, Year)  %>%
  unique()%>%
  arrange(PLOT, SUBP, Year) |>
  mutate(year0 = lag(Year), .by = c(PLOT, SUBP)) |>
  drop_na(year0) 

# need to get rid of the older sites adn the sites with disturban not including

# ---- 2. Tree-level components per plot-interval -------------------------------
live <- tree.demo.subp |> filter(status == "live") |> select(PLOT, SUBP, Year, UniqueID, above_dbh_biomass)

t1 <- intervals |> inner_join(live |> rename(year0 = Year, b1 = above_dbh_biomass), by = c("PLOT", "SUBP","year0"))
t2 <- intervals |> inner_join(live |> rename(b2 = above_dbh_biomass),               by = c("PLOT", "SUBP","Year"))


comp <- full_join(t1, t2, by = c("PLOT", "SUBP","year0", "Year", "UniqueID")) |>
  summarise(
    growth  = sum(b2 - b1, na.rm = TRUE),          # survivors only (NA otherwise)
    recruit = sum(b2[is.na(b1)]),                  # ingrowth (incl. its growth)
    mort    = sum(b1[is.na(b2)]),                  # biomass at t1 of trees that died
    net     = sum(b2, na.rm = TRUE) - sum(b1, na.rm = TRUE),
    dbh_death =mean(b1[is.na(b2)]),
    mean_dbh = mean(b2, na.rm = T),
    .by = c(PLOT, SUBP, year0, Year)
  ) |>
  right_join(intervals, by = c("PLOT", "SUBP", "year0", "Year")) |>
  mutate(across(c(growth, recruit, mort, net), \(x) replace_na(x, 0)/100)) %>%
  mutate(dbh_death = replace_na(dbh_death, 0)) %>%
  left_join(subplot_traj) %>%
  mutate(
    interval  = Year - year0,
    across(c(growth, recruit, mort, net), \(x) x / interval, .names = "{.col}_ann")) %>%
  left_join(ages) %>%
  mutate(age2 = Year - age)

ggplot(comp, aes(recruit_ann, mort_ann, color =  log(age2)))+
  geom_point() +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

ggplot(comp, aes(net_ann, mort_ann, color = decid.index.msq.sp))+
#ggplot(comp, aes(net_ann, mort_ann, color = log(age2)))+ 
  geom_point() +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

ggplot(comp, aes(net, recruit, color = log(age2)))+
  geom_point() +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

ggplot(comp, aes(net_ann, growth_ann, color = log(age2)))+
  geom_point() +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

ggplot(comp, aes(growth_ann, mort_ann, color = log(age2)))+
  geom_point() +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())


# subset to younger
young <- comp %>%
  filter(age2 < 75)
ggplot(young, aes(decid.index.msq.sp, mort_ann, color = age2))+
  geom_point() +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())


ggplot(young, aes(mean_dbh, dbh_death, color = decid.index.msq.sp))+
  geom_point() +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40")



# Prop of dead deciduous and density of live spruce
live <- tree.demo.subp |> filter(status == "live") |> select(PLOT, SUBP, Year, SPECIES, func_gp, UniqueID, above_dbh_biomass)
t1 <- intervals |> inner_join(live |> rename(year0 = Year, b1 = above_dbh_biomass), by = c("PLOT", "SUBP","year0"))
t2 <- intervals |> inner_join(live |> rename(b2 = above_dbh_biomass),               by = c("PLOT", "SUBP","Year"))

comp2 <- full_join(t1, t2, by = c("PLOT", "SUBP","year0", "Year", "UniqueID", "func_gp")) |>
  summarise(
    growth  = sum(b2 - b1, na.rm = TRUE),          # survivors only (NA otherwise)
    recruit = sum(b2[is.na(b1)]),                  # ingrowth (incl. its growth)
    mort    = sum(b1[is.na(b2)]),                  # biomass at t1 of trees that died
    net     = sum(b2, na.rm = TRUE) - sum(b1, na.rm = TRUE),
    dbh_death =mean(b1[is.na(b2)]),
    mean_dbh = mean(b2, na.rm = T),
    .by = c(PLOT, SUBP, year0, Year, func_gp)
  ) |>
  right_join(intervals, by = c("PLOT", "SUBP", "year0", "Year")) |>
  mutate(across(c(growth, recruit, mort, net), \(x) replace_na(x, 0)/100)) %>%
  mutate(dbh_death = replace_na(dbh_death, 0)) %>%
  left_join(subplot_traj) %>%
  mutate(
    interval  = Year - year0,
    across(c(growth, recruit, mort, net), \(x) x / interval, .names = "{.col}_ann")) %>%
  left_join(ages) %>%
  mutate(age2 = Year - age)

comp_wide <- comp2 |>
  select(-c(growth, recruit, mort, net)) %>%
  filter(!is.na(func_gp)) |>                                   # drop empty-subplot rows from right_join
  pivot_wider(id_cols     = c(PLOT, SUBP, year0, Year),
              names_from  = func_gp,
              values_from = c(growth_ann, recruit_ann, mort_ann, net_ann, dbh_death, mean_dbh),
              names_glue  = "{func_gp}_{.value}",
              values_fill = 0) |>
  right_join(intervals, by = c("PLOT", "SUBP", "year0", "Year")) 

glimpse(comp_wide)

ggplot(comp_wide, aes(Decid_Broad_mort_ann, Conifer_growth_ann, color = Decid_Broad_growth_ann))+
  geom_point()

# select sites where there is high death
glimpse(comp)

death <- comp %>%
  filter(mort_ann>0) %>%
  filter(decid.index.msq.sp > 0)

comp.wide.death <- comp_wide %>%
  filter(PLOT %in% death$PLOT) %>%
  # left_join(tree.summary) %>%
  # left_join(tree.summary.species, by = c("PLOT", "Year", "area.sampled")) %>%
  left_join(comp)%>%
  mutate(plot_subp = paste(PLOT, SUBP, sep = "_")) %>%
  left_join(ages) %>%
  mutate(age2 = Year - age) %>%
  left_join(tree.subplot)

glimpse(comp.wide.death) 

ggplot(comp.wide.death, aes(Decid_Broad_dens, Decid_Broad_mort_ann, color = log(age2))) +
  geom_point()+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())


ggplot(comp.wide.death, aes(age2, Decid_Broad_mort_ann, color = agb_m2)) +
  geom_point()+
  xlim(10, 75)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())


subplot_sp <- trees.2 %>%
  mutate(above_dbh_biomass = allom.dbh(SPECIES, DBH),
         BA = 0.00007854 * (DBH)^2) %>%
  mutate(func_gp = ifelse(SPECIES %in% decid, "Decid_Broad", "Conifer")) %>%
  filter(PLOT %!in% c("BDI5", "BDI6", "FP5C", "FP5A", "FP5D")) %>%
  filter(Network == "RSN") %>%
  filter(DA == 1) %>%
  group_by(Year, PLOT, SUBP, func_gp) %>%
  summarize(biomass = sum(above_dbh_biomass), ba = sum(BA), n=n()) %>%
  ungroup()%>%
  mutate(biomass_m2_sp = biomass/100, ba_m2_sp = ba/100, dens_sp = n/100) %>%
  select(-c(biomass, ba, n)) %>%
  pivot_wider(id_cols = c(Year, PLOT, SUBP),
              names_from = func_gp,
              values_from = c(dens_sp, ba_m2_sp, biomass_m2_sp),
              names_glue = "{func_gp}_{.value}") %>%
  mutate(across(ends_with(c("_dens_sp", "_ba_m2_sp", "biomass_m2_sp")), ~replace_na(.x, 0))) %>% 
  mutate(rel_dens_decid = Decid_Broad_dens_sp/(Decid_Broad_dens_sp + Conifer_dens_sp), 
         rel_biomass_decid = Decid_Broad_biomass_m2_sp/(Decid_Broad_biomass_m2_sp + Conifer_biomass_m2_sp), 
         decid.index.msq.sp = (rel_dens_decid + rel_biomass_decid)/2) %>%
  select(-c(rel_dens_decid, rel_biomass_decid)) %>%
  mutate(decid.index.msq.sp = replace_na(decid.index.msq.sp, 0))

comp.wide.death2 <- left_join(comp.wide.death, subplot_sp) %>%
  filter(age2 > 40) %>%
  filter(Decid_Broad_mort_ann > 0)

summary(comp.wide.death2$age2)
glimpse(comp.wide.death2)

ggplot(comp.wide.death2, aes(Conifer_dens_sp, Decid_Broad_mort_ann, color = Decid_Broad_dens_sp)) +
  geom_point()+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())


ggplot(comp.wide.death2, aes(Decid_Broad_dens_sp, Conifer_mort_ann, color = Decid_Broad_dens_sp)) +
  geom_point()+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

ggplot(comp.wide.death2, aes(dens_sp, mort_ann, color = (Decid_Broad_dens_sp/dens_sp))) +
  geom_point()+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

ggplot(comp.wide.death2, aes(dens_sp, mort_ann, color = (Decid_Broad_dens_sp/dens_sp))) +
  geom_point()+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

ggplot(comp.wide.death2, aes(dens_sp, mort_ann, color = decid.index.msq.sp)) +
  geom_point()+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

ggplot(comp.wide.death2, aes(dens_sp, mort_ann, color = dbh_death)) +
  geom_point()+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())


## By the subplot - i want numbers of trees that are dying and by species
live <- tree.demo.subp |> filter(status == "live") |> select(PLOT, SUBP, Year, SPECIES, func_gp, UniqueID, above_dbh_biomass)
t1 <- intervals |> inner_join(live |> rename(year0 = Year, b1 = above_dbh_biomass), by = c("PLOT", "SUBP","year0"))
t2 <- intervals |> inner_join(live |> rename(b2 = above_dbh_biomass),               by = c("PLOT", "SUBP","Year"))

comp.n <- full_join(t1, t2, by = c("PLOT", "SUBP", "year0", "Year", "UniqueID", "func_gp")) |>
  summarise(
    growth     = sum(b2 - b1, na.rm = TRUE),
    recruit    = sum(b2[is.na(b1)]),
    mort       = sum(b1[is.na(b2)]),
    net        = sum(b2, na.rm = TRUE) - sum(b1, na.rm = TRUE),
    dbh_death  = mean(b1[is.na(b2)]),
    mean_dbh   = mean(b2, na.rm = TRUE),
    n_dead     = sum(!is.na(b1) &  is.na(b2)),
    n_recruit  = sum( is.na(b1) & !is.na(b2)),
    n_survivor = sum(!is.na(b1) & !is.na(b2)),
    n_t1       = sum(!is.na(b1)),
    n_t2       = sum(!is.na(b2)),
    .by = c(PLOT, SUBP, year0, Year, func_gp)
  ) |>
  pivot_wider(names_from  = func_gp,
              values_from = growth:n_t2,                     # every summarized column
              names_glue  = "{func_gp}_{.value}",
              values_fill = 0) |>
  right_join(intervals, by = c("PLOT", "SUBP", "year0", "Year")) |>   # empty subplots back in
  mutate(across(where(is.numeric), \(x) replace_na(x, 0))) |>
  left_join(subplot_traj, by = c("PLOT", "SUBP", "Year")) |>
  mutate(
    interval = Year - year0,
    across(matches("_(growth|recruit|mort|net)$"), \(x) x / interval, .names = "{.col}_ann")
  ) |>
  left_join(ages) |>
  mutate(age2 = Year - age)

# reduce the dataset to decid intermediate
glimpse(comp.n)
comp.n2 <- comp.n %>%
  filter(age2 < 70) %>%
  filter(age2 > 35) %>%
  filter(decid.index.msq.sp>0) %>%
  filter(PLOT != "MDI1")

glimpse(comp.n2)  


ggplot(comp.n, aes(Decid_Broad_n_dead, Conifer_n_t2 - Conifer_n_t1))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())


ggplot(comp.n2, aes(Decid_Broad_n_dead, Conifer_n_t2 + Decid_Broad_n_t2))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

ggplot(comp.n2, aes(Decid_Broad_n_t2, Conifer_n_t2 - Conifer_n_t1))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())


ggplot(comp.n2, aes(Decid_Broad_n_dead, Conifer_n_t2))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())


ggplot(comp.n2, aes(Decid_Broad_n_t2, Conifer_n_t2))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  ylim(0, 150)+
  xlim(0, 150)

ggplot(comp.n2, aes(Decid_Broad_n_t1, Conifer_n_t1))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())+
  ylim(0, 150)+
  xlim(0, 150)

ggplot(comp.n, aes(Decid_Broad_n_t2 - Decid_Broad_n_t1, Conifer_n_t2-Conifer_n_t1, color = decid.index.msq.sp))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

ggplot(comp.n2, aes(Decid_Broad_n_t2 - Decid_Broad_n_t1, Conifer_n_t2-Conifer_n_t1, color = age2))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())


ggplot(comp.n2, aes(Decid_Broad_n_dead, Decid_Broad_n_t1))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

ggplot(comp.n2, aes(Decid_Broad_n_dead, Decid_Broad_n_t2 - Decid_Broad_n_t1))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

ggplot(comp.n2, aes(Decid_Broad_n_dead, Decid_Broad_n_t2 + Conifer_n_t2))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

ggplot(comp.n2, aes(Conifer_n_dead, Decid_Broad_n_t2 - Decid_Broad_n_t1))+
  geom_point() +
  theme_bw(base_size = 14) +
  geom_hline(yintercept = 0)+
  theme(panel.grid = element_blank())


# comp n age categorie
comp.n3 <- comp.n %>%
  mutate(age_class = case_when(
    age2 < 40 ~ "Young", 
    age2 > 80 ~ "Mature",
    TRUE ~ "Intermediate"
  )) %>%
  filter(age_class == "Intermediate")

glimpse(comp.n3)
ggplot(comp.n3, aes(Decid_Broad_n_t2 - Decid_Broad_n_t1, Conifer_n_t2-Conifer_n_t1, color = decid.index.msq.sp))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~age_class, scales = "free")


ggplot(comp.n3, aes(Decid_Broad_net_ann, Conifer_net_ann, color = decid.index.msq.sp))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~age_class, scales = "free")

ggplot(comp.n3, aes(Decid_Broad_net_ann, Conifer_net_ann, color = Decid_Broad_n_t2 + Conifer_n_t2))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~age_class, scales = "free")

ggplot(comp.n3, aes(Decid_Broad_net_ann, Conifer_net_ann, color = Decid_Broad_n_t1))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~age_class, scales = "free")

ggplot(comp.n3, aes(Decid_Broad_net_ann, Conifer_net_ann, color = Conifer_n_t1))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~age_class, scales = "free")

ggplot(comp.n3, aes(Decid_Broad_n_t2 - Decid_Broad_n_t1, Conifer_n_t2-Conifer_n_t1, color = log(Decid_Broad_dbh_death)))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~age_class, scales = "free")

ggplot(comp.n3, aes(Decid_Broad_n_t2 - Decid_Broad_n_t1, Conifer_n_t2-Conifer_n_t1, color = age_class))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 0)
  #facet_wrap(~age_class, scales = "free")

ggplot(comp.n3, aes(Decid_Broad_n_t2 - Decid_Broad_n_t1, Conifer_n_t2-Conifer_n_t1, color = Decid_Broad_n_t1))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~age_class, scales = "free")

ggplot(comp.n3, aes(Decid_Broad_n_t2 - Decid_Broad_n_t1, Conifer_n_t2-Conifer_n_t1, color = Decid_Broad_n_t2))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~age_class, scales = "free")

ggplot(comp.n3, aes(Decid_Broad_n_t2 - Decid_Broad_n_t1, Conifer_n_t2-Conifer_n_t1, color = Decid_Broad_n_t2 + Conifer_n_t2))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~age_class, scales = "free")

ggplot(comp.n3, aes(age2, Conifer_recruit_ann, color = decid.index.msq.sp))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~age_class, scales = "free")

ggplot(comp.n3, aes(age2, Decid_Broad_recruit_ann, color = decid.index.msq.sp))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~age_class, scales = "free")
ggplot(comp.n3, aes(age2, Conifer_n_recruit_ann, color = decid.index.msq.sp))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())
Conifer_n_recruit_ann

ggplot(comp.n3, aes(age2, Decid_Broad_n_dead+Conifer_n_dead, color = decid.index.msq.sp))+
  geom_point() +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

glimpse(comp.n3)
# sites where number of confier trees are decreasing
dec.con <- comp.n3 %>%
  filter(Conifer_n_t2-Conifer_n_t1 < 0) %>%
  group_by(PLOT) %>%
  summarize(n = length(unique(SUBP)))

unique(dec.con$PLOT)

check <- comp.n3 %>%
  mutate(Conifer_n_t2-Conifer_n_t1 )

# age and no of trees
glimpse(tree.combo)
ggplot(tree.combo, aes(age2, Conifer_dens, color = decid.index.msq.sp))+
  geom_point(alpha = 0.2) +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  xlim(0, 75)+
  theme(panel.grid = element_blank()) 

ggplot(tree.combo, aes(age2, Decid_Broad_dens, color = decid.index.msq.sp))+
  geom_point(alpha = 0.2) +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  xlim(0, 75)+
  theme(panel.grid = element_blank()) 


ggplot(tree.combo, aes(age2, Conifer_dens, color = decid.index.msq.sp))+
  geom_point(alpha = 0.2) +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  xlim(0, 75)+
  theme(panel.grid = element_blank()) 

ggplot(tree.combo, aes(age2, Decid_Broad_dens, color = decid.index.msq.sp))+
  geom_point(alpha = 0.2) +
  geom_hline(yintercept = 0)+
  theme_bw(base_size = 14) +
  xlim(0, 75)+
  theme(panel.grid = element_blank()) 

# Model comp
comp.n4 <- comp.n %>%
  mutate(age_class = case_when(
    age2 < 40 ~ "Young", 
    age2 > 80 ~ "Mature",
    TRUE ~ "Intermediate"
  )) %>%
  filter(age_class != "Mature") %>%
  left_join(pc1)




# Conifer rec vs. age and decid index

# Conifer biomass over time


library(tidyverse)
library(glmmTMB)
library(DHARMa)

# ---- Data prep: one row per subplot, covariates at start of interval ----------
sp <- comp.n4 |>
  mutate(PLOT     = factor(PLOT),
         age0     = age2 - interval,                     # age at first census (check: age2 = end age)
         age_c    = (age0 - mean(age0)) / sd(age0),      # scaled for convergence
         decid0   = decid.index.msq.sp,                  # subplot composition at t
         topo_c = topo_moist_PC1,
         area_m2  = 100)

# ---- 1. Conifer recruitment vs. age x deciduous canopy (relay test) ----------
# Counts of new conifer stems; offset makes it a rate per m2 per year
m_rec_n <- glmmTMB(Conifer_n_recruit ~ age_c * decid0 + topo_c +
                     offset(log(area_m2 * interval)) + (1 | PLOT),
                   family = nbinom2, data = sp)
# Relay prediction: positive age_c:decid0 (more conifer recruitment under older deciduous canopies)

# Biomass version (zeros + skew)
m_rec_b <- glmmTMB(Conifer_recruit_ann ~ age_c * decid0 + topo_c + (1 | PLOT),
                   family = tweedie(link = "log"), data = sp)

# ---- 2. Where do conifer stems come from: recruits vs survivors -------------
# Proportion of conifer stems alive at t2 that are new recruits
m_share <- glmmTMB(cbind(Conifer_n_recruit, Conifer_n_survivor) ~ age_c * decid0 + topo_c +
                     (1 | PLOT),
                   family = betabinomial, data = filter(sp, Conifer_n_t2 > 0))
# Relay: share stays high / rises with age. Initial floristics: share declines toward 0.

# ---- 3. Deciduous mortality hazard vs. age ---------------------------------
# Proportion of deciduous stems alive at t1 that died; cloglog + log(interval)
# offset gives an annual hazard, so differing intervals are handled properly
m_dmort <- glmmTMB(cbind(Decid_Broad_n_dead, Decid_Broad_n_survivor) ~ age_c * decid0 + topo_c +
                     offset(log(interval)) + (1 | PLOT),
                   family = betabinomial(link = "cloglog"), data = filter(sp, Decid_Broad_n_t1 > 0))

# ---- Checks ----------------------------------------------------------------
mods <- list(rec_n = m_rec_n, rec_b = m_rec_b, share = m_share, dmort = m_dmort)
map(mods, summary)
iwalk(mods, \(m, y) plot(simulateResiduals(m), title = y))
testZeroInflation(simulateResiduals(m_rec_n))   # if significant: add ziformula = ~ 1

library(performance)
r2(m_rec_n)                # Nakagawa R2 (marginal + conditional)
icc(m_rec_n)               # share of non-fixed variance due to plots
m_q   <- update(m_rec_n, . ~ . - age_c * decid0 + poly(age_c, 2) * decid0)
m_log <- update(m_rec_n, . ~ . - age_c * decid0 + log(age0) * decid0)
AIC(m_rec_n, m_q, m_log)

res_q <- simulateResiduals(m_q)
plot(res_q)
plotResiduals(res_q, sp$age0, main = "age")
testDispersion(res_q)



library(marginaleffects)

age_mu <- mean(sp$age0); age_sd <- sd(sp$age0)

nd <- expand_grid(age0 = seq(min(sp$age0), max(sp$age0), length.out = 100),
                  decid0 = c(0, .5, 1)) |>
  mutate(age_c = (age0 - age_mu) / age_sd,
         topo_c = 0, area_m2 = 10000, interval = 1, PLOT = NA)   # offset -> per ha per yr

pred <- predictions(m_rec_n, newdata = nd, re.form = NA) |> as_tibble()

decid_pal <- setNames(colorRampPalette(c("#1b4d3e", "#7fa37a", "#e8d9a0", "#d4a017", "#b5651d"))(5), c("0","0.25", "0.5", "0.75", "1"))

ggplot(pred, aes(age0, estimate, colour = factor(decid0), fill = factor(decid0))) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_rug(data = sp, aes(x = age0), inherit.aes = FALSE, alpha = .2) +
  scale_colour_manual(values = decid_pal, aesthetics = c("colour", "fill"),
                      labels = c("0%", "50%", "100%")) +
  labs(x = "Stand age at first census (yr)",
       y = expression(Conifer~recruitment~(stems~ha^-1~yr^-1)),
       colour = "Deciduous", fill = "Deciduous") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) 
