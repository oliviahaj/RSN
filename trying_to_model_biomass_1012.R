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
  select(c(Network, Year, PLOT, decid.index.msq, sample.pd, age2))%>%
  filter(Network == "RSN")

# look at decid index and age
ggplot(tree.sum2, aes(age2, decid.index.msq, color = Network))+
  geom_point()+
  theme_bw()

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
glimpse(tree.miss)

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

# tree and biomass
tree.biomass <- trees.2 %>%
  mutate(above_dbh_biomass = allom.dbh(SPECIES, DBH),
         live = ifelse(DA == 0, "dead", "live")) %>%
  filter(Network == "RSN") %>%
  filter(PLOT %!in% c("BDI5", "BDI6", "FP5C", "FP5A", "FP5D"))

# read in site productivity and join in
drivers <- read.csv('/Users/olhajek/Desktop/RSN/RSN_proj/Data/2026 Data/RSN_Site/RSN_Climate_Site_characteristcs.csv') 

pc1 <- drivers %>%
  select(c(PLOT, topo_moist_PC1)) %>%
  unique()

tree.meta2 <- left_join(tree.meta, pc1) %>%
  filter(age2 < 80) %>%
  filter(Network == "RSN") %>%
  filter(PLOT %!in% c("BDI5", "BDI6", "FP5C", "FP5A", "FP5D")) %>%
  distinct()
  

# join everything
treess <- left_join(tree.biomass, pc1)
trees.1 <- left_join(treess, tree.meta)
glimpse(trees.1)
tree.select <- trees.1 %>%
  select(c(UniqueID, Year, PLOT, DA, live, above_dbh_biomass, decid.index.msq, sample.pd, age2))%>%
  rename(status = live) %>%
  filter(age2 < 80)

# We're going to subset to RSN and 70 years....i think
# could increase wiht CAFI, but issues iwth diameter, white spruce vs. black spruce...index might need to shift?
# maybe that gets at teh site quality with the PC??

## models
set.seed(42)
decid_levels <- c(0, .25, .5, .75, 1)
n_draws      <- 1000

# ---- 1. Census intervals per plot ---------------------------------------------
intervals <- tree.meta2 |>
  arrange(PLOT, Year) |>
  mutate(year0 = lag(Year), .by = PLOT) |>
  drop_na(year0) |>
  select(PLOT, year0, Year) 

# need to get rid of the older sites adn the sites with disturban not including

# ---- 2. Tree-level components per plot-interval -------------------------------
live <- tree.select |> filter(status == "live") |> select(PLOT, Year, UniqueID, above_dbh_biomass)

t1 <- intervals |> inner_join(live |> rename(year0 = Year, b1 = above_dbh_biomass), by = c("PLOT", "year0"))
t2 <- intervals |> inner_join(live |> rename(b2 = above_dbh_biomass),               by = c("PLOT", "Year"))

t1 |> count(PLOT, year0, Year, UniqueID) |> filter(n > 1)
t2 |> count(PLOT, year0, Year, UniqueID) |> filter(n > 1)


plot_area_m2 <- 1200          # 12 x 100 m2 subplots
kg_to_Mg_ha  <- (10000 / plot_area_m2) / 1000   # kg per plot -> Mg/ha

comp <- full_join(t1, t2, by = c("PLOT", "year0", "Year", "UniqueID")) |>
  summarise(
    growth  = sum(b2 - b1, na.rm = TRUE),          # survivors only (NA otherwise)
    recruit = sum(b2[is.na(b1)]),                  # ingrowth (incl. its growth)
    mort    = sum(b1[is.na(b2)]),                  # biomass at t1 of trees that died
    net     = sum(b2, na.rm = TRUE) - sum(b1, na.rm = TRUE),
    .by = c(PLOT, year0, Year)
  ) |>
  right_join(intervals, by = c("PLOT", "year0", "Year")) |>
  mutate(across(c(growth, recruit, mort, net), \(x) replace_na(x, 0) * kg_to_Mg_ha)) 

# ---- 3. Annualize + start-of-interval covariates ------------------------------
glimpse(comp)
glimpse(tree.meta2)
tree.meta2 |> count(PLOT, Year) |> filter(n > 1)
comp |> count(PLOT, year0) |> filter(n > 1)
comp |> count(PLOT, year0, Year) |> filter(n > 1) 

dat <- comp |>
  left_join(tree.meta2 |> select(PLOT, year0 = Year, age2, decid.index.msq, topo_moist_PC1),
            by = c("PLOT", "year0")) |>
  mutate(
    interval  = Year - year0,
    across(c(growth, recruit, mort, net), \(x) x / interval, .names = "{.col}_ann"),
    age_mid   = age2 + interval / 2,
    log_int_c = log(interval) - mean(log(interval)),   # census-interval bias term
    site_q_c  = topo_moist_PC1,   # 0 = average site quality
    PLOT = factor(PLOT)
  )

# ---- 4. Data checks -----------------------------------------------------------
dat |> ggplot(aes(age_mid, decid.index.msq)) + geom_point() + theme_bw()   # coverage / confounding
cor(dat$age_mid, dat$decid.index.msq)
dat |> summarise(across(ends_with("_ann"), list(zero = \(x) mean(x == 0), neg = \(x) mean(x < 0))))

# ---- 5. Model specs -----------------------------------------------------------
# Growth: Gamma if strictly positive; otherwise Gaussian (negative growth from breakage etc.)
fams <- list(
  growth_ann  = if (all(dat$growth_ann > 0)) Gamma(link = "log") else gaussian(),
  recruit_ann = tweedie(link = "log"),   # many zeros, non-negative
  mort_ann    = tweedie(link = "log"),
  net_ann     = gaussian()               # can be negative
)

rhs <- "ns(age_mid, df = 3) * decid.index.msq + log_int_c + (1 | PLOT)"
fit_all <- function(rhs, data) imap(fams, \(fam, y)
                                    glmmTMB(as.formula(paste(y, "~", rhs)), data = data, family = fam))

# ---- 6. Does site quality earn its place? (compare on the same rows) ----------
dat_sq   <- dat |> drop_na(topo_moist_PC1)
mods_sq  <- fit_all(paste(rhs, "+ topo_moist_PC1"), dat_sq)
mods_nsq <- fit_all(rhs, dat_sq)
map2_dfr(mods_sq, mods_nsq, \(a, b) tibble(dAIC = AIC(a) - AIC(b)), .id = "component")
map(mods_sq, \(m) summary(m)$coefficients$cond["topo_moist_PC1", ])

use_sq <- TRUE   # <- set after reviewing the above
mods <- if (use_sq) mods_sq else fit_all(rhs, dat)   # refit on all rows if dropping site_q

# ---- 7. Diagnostics -----------------------------------------------------------
library(DHARMa)
walk(mods, \(m) plot(simulateResiduals(m)))
map(mods, summary)

res <- map(mods, simulateResiduals)

iwalk(mods, \(m, y) {
  mf <- model.frame(m)                       # rows actually used in the fit
  par(mfrow = c(2, 2))
  plotResiduals(res[[y]], mf$age_mid,     main = paste(y, "- age"))
  plotResiduals(res[[y]], mf$decid_index, main = paste(y, "- decid"))
  plotResiduals(res[[y]], mf$log_int_c,   main = paste(y, "- interval"))
  plotResiduals(res[[y]], mf$site_q_c,    main = paste(y, "- moisture"))
  par(mfrow = c(1, 1))
})

iwalk(res, \(r, y) { cat("\n", y, "\n"); print(testZeroInflation(r)); print(testDispersion(r)) })

# ---- 8. Prediction grid: annual ages, fixed composition, avg site & interval --
a0   <- floor(min(dat$age_mid))     # set to 0 to anchor at the fire (extrapolates)
amax <- ceiling(max(dat$age_mid))

nd <- expand_grid(decid.index.msq = decid_levels, age_mid = seq(a0, amax, by = 1)) |>
  mutate(topo_moist_PC1 = 0, log_int_c = 0)

# Starting stock: mean observed biomass in the youngest stands (or set B0 <- 0)
B0 <- tree.select |> filter(age2 <= a0 + 5) |> summarise(mean(above_dbh_biomass, na.rm = TRUE)) |> pull()

# ---- 9. Simulate rate curves from fixed-effect uncertainty --------------------
sim_rates <- function(m, nd, n = n_draws) {
  X <- model.matrix(delete.response(terms(m)), nd)            # reuses fitted ns() knots
  B <- MASS::mvrnorm(n, fixef(m)$cond, vcov(m)$cond)
  family(m)$linkinv(X %*% t(B))                               # rows = grid, cols = draws
}

long <- imap(mods, \(m, y)
             nd |> select(decid.index.msq, age_mid) |>
               bind_cols(as_tibble(sim_rates(m, nd), .name_repair = \(x) paste0("draw_", seq_along(x)))) |>
               pivot_longer(starts_with("draw_"), names_to = "draw", values_to = y)) |>
  reduce(left_join, by = c("decid.index.msq", "age_mid", "draw"))

# ---- 10. Integrate rates over age (within each draw) --------------------------
cum <- long |>
  arrange(decid.index.msq, draw, age_mid) |>
  mutate(
    cum_production = cumsum(growth_ann + recruit_ann),
    cum_mortality  = cumsum(mort_ann),
    stock_comp     = B0 + cum_production - cum_mortality,   # stock from components
    stock_net      = B0 + cumsum(net_ann),                  # stock from net model
    .by = c(decid.index.msq, draw)
  )

# ---- 11. Summarise draws (rates + stocks) -------------------------------------
summ <- cum |>
  pivot_longer(-c(decid.index.msq, age_mid, draw), names_to = "metric") |>
  summarise(est = median(value),
            lo  = quantile(value, .025),
            hi  = quantile(value, .975),
            .by = c(decid.index.msq, age_mid, metric))

# ---- 12. Truncate each composition level to its observed age range ------------
support <- map_dfr(decid_levels, \(d)
                   dat |> filter(abs(decid.index.msq - d) <= .125) |>
                     summarise(decid.index.msq = d, age_lo = min(age_mid), age_hi = max(age_mid)))

summ_obs <- summ |>
  left_join(support, by = "decid.index.msq") |>
  filter(age_mid <= age_hi)            # keep from a0 so integration start is shared

# ---- 13. Plots ----------------------------------------------------------------
pal <- list(scale_color_viridis_d(end = .9, aesthetics = c("color", "fill")), theme_bw())
lab_decid <- \(x) factor(scales::percent(x), levels = scales::percent(decid_levels))

# (a) Annual rates
summ_obs |>
  filter(str_ends(metric, "_ann")) |>
  mutate(decid = lab_decid(decid.index.msq),
         metric = fct_recode(metric, Growth = "growth_ann", Recruitment = "recruit_ann",
                             Mortality = "mort_ann", Net = "net_ann")) |>
  ggplot(aes(age_mid, est, color = decid, fill = decid)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .15, color = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = "Stand age (yr)", y = expression(Biomass~flux~(Mg~ha^-1~yr^-1)),
       color = "Deciduous", fill = "Deciduous") + pal

# (b) Integrated stock curves with observed plot trajectories
summ_obs |>
  filter(metric == "stock_comp") |>                        # swap to "stock_net" to compare
  mutate(decid = lab_decid(decid.index.msq)) |>
  ggplot(aes(age_mid, est, color = decid, fill = decid)) +
  # geom_line(data = tree.select, aes(age2, above_dbh_biomass, group = PLOT),
  #           inherit.aes = FALSE, color = "grey75", alpha = .6) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .15, color = NA) +
  geom_line(linewidth = 1) +
  labs(x = "Stand age (yr)", y = expression(Aboveground~biomass~(Mg~ha^-1)),
       color = "Deciduous", fill = "Deciduous") + pal

# (c) Cumulative production vs. mortality
summ_obs |>
  filter(metric %in% c("cum_production", "cum_mortality")) |>
  mutate(decid = lab_decid(decid.index.msq),
         metric = fct_recode(metric, Production = "cum_production", Mortality = "cum_mortality")) |>
  ggplot(aes(age_mid, est, color = decid, linetype = metric)) +
  geom_line(linewidth = 1) +
  labs(x = "Stand age (yr)", y = expression(Cumulative~biomass~(Mg~ha^-1)),
       color = "Deciduous", linetype = NULL) + pal




library(tidyverse)
library(mgcv)

decid_levels <- c(0, .25, .5, .75, 1)
n_draws      <- 1000
set.seed(42)

# ---- 1. Long data: one row per plot-interval x process ----------------------
procs <- dat |>
  drop_na(age_mid, decid.index.msq, log_int_c, topo_moist_PC1) |>
  mutate(PLOT       = factor(PLOT),
         topo_c     = topo_moist_PC1 - mean(topo_moist_PC1),   # 0 = average moisture
         growth_ann = pmax(growth_ann, 0)) |>                  # Tweedie needs >= 0
  select(PLOT, age_mid, decid.index.msq, log_int_c, topo_c,
         Growth = growth_ann, Recruit = recruit_ann, Death = mort_ann, Net = net_ann) |>
  pivot_longer(c(Growth, Recruit, Death, Net), names_to = "demo_process", values_to = "rate")

# ---- 2. Common age range supported at every composition level ---------------
rng <- map_dfr(decid_levels, \(dl)
               procs |> distinct(PLOT, age_mid, decid.index.msq) |>
                 filter(abs(decid.index.msq - dl) <= .125) |>
                 summarise(decid = dl, n = n(), lo = min(age_mid), hi = max(age_mid)))
rng                                                # check n > 0 for every level
ages <- seq(ceiling(max(rng$lo)), floor(min(rng$hi)), by = 1)

# ---- 3. Model + integration functions ---------------------------------------
fams <- list(Growth = tw(), Recruit = tw(), Death = tw(), Net = gaussian())

fit_proc <- function(d, proc) {
  gam(rate ~ s(age_mid, k = 5) + s(age_mid, by = decid.index.msq, k = 5) +
        log_int_c + topo_c + s(PLOT, bs = "re"),
      data = d, family = fams[[proc]], method = "REML")
}

cum_draws <- function(m, proc, plot_ref, ages, n_draws = 1000) {
  nd <- expand_grid(decid = decid_levels, age_mid = ages) |>
    mutate(decid.index.msq = decid, log_int_c = 0, topo_c = 0, PLOT = plot_ref)
  X    <- predict(m, nd, type = "lpmatrix", exclude = "s(PLOT)")   # population-level
  beta <- MASS::mvrnorm(n_draws, coef(m), vcov(m, unconditional = TRUE))
  sims <- m$family$linkinv(X %*% t(beta))                          # response scale
  if (proc == "Death") sims <- -sims                               # losses plot below 0
  
  nd |> select(decid, age_mid) |>
    bind_cols(as_tibble(sims, .name_repair = \(x) str_c("draw_", seq_along(x)))) |>
    pivot_longer(starts_with("draw_"), names_to = "draw", values_to = "rate") |>
    arrange(decid, draw, age_mid) |>
    mutate(cum = cumsum(c(0, diff(age_mid) * (head(rate, -1) + tail(rate, -1)) / 2)),
           .by = c(decid, draw)) |>                                # trapezoid integral
    summarise(rate = mean(rate), est = mean(cum),
              lo = quantile(cum, .025), hi = quantile(cum, .975),
              .by = c(decid, age_mid))
}

# ---- 4. Fit, check, integrate -----------------------------------------------
fits <- procs |>
  nest(.by = demo_process) |>
  mutate(model = map2(data, demo_process, fit_proc),
         cum   = map2(model, demo_process, \(m, p) cum_draws(m, p, procs$PLOT[1], ages, n_draws)))

walk2(fits$model, fits$demo_process, \(m, p) { cat("\n====", p, "====\n"); print(summary(m)) })
# walk(fits$model, gam.check)   # basis size + convergence

cum_proc <- fits |> select(demo_process, cum) |> unnest(cum) |>
  mutate(decid_lab = factor(str_c(decid * 100, "% deciduous"),
                            levels = str_c(decid_levels * 100, "% deciduous")))

# ---- 5. Budget figure -------------------------------------------------------
proc_col <- c(Growth = "#1b9e77", Recruit = "#7570b3", Death = "#d95f02")

comp  <- cum_proc |> filter(demo_process != "Net") |>
  mutate(demo_process = factor(demo_process, levels = c("Recruit", "Growth", "Death")))
net_l <- cum_proc |> filter(demo_process == "Net")
sum_l <- comp |> summarise(est = sum(est), .by = c(decid_lab, age_mid))

p_budget <- ggplot() +
  geom_area(data = comp, aes(age_mid, est, fill = demo_process), alpha = 0.75) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_ribbon(data = net_l, aes(age_mid, ymin = lo, ymax = hi), fill = "black", alpha = 0.15) +
  geom_line(data = sum_l, aes(age_mid, est), linetype = "dashed", colour = "grey20") +
  geom_line(data = net_l, aes(age_mid, est), linewidth = 1.4, colour = "black") +
  facet_wrap(~ decid_lab, nrow = 1) +
  scale_fill_manual(values = proc_col) +
  labs(x = "Stand age (years)",
       y = expression(Cumulative~AGB~flux~(Mg~ha^-1)),
       fill = NULL,
       caption = "Solid line ± 95% CI = Net; dashed = Growth + Recruit + Death") +
  theme_bw() +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_text(face = "bold"))
p_budget

