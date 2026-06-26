subs <- read.csv("/Users/olhajek/Desktop/RSN/CAFI_data/452_CAFI_SUBPLOT_v20241220.csv")

library(readxl)
mel_age <- read_excel("/Users/olhajek/Desktop/RSN/CAFI_data/CAFI_Boyd_061522/CAFI_stand_age/age_cafi_021621.xlsx", sheet = 1)
trees <- read_excel('/Users/olhajek/Desktop/RSN/CAFI_data/CAFI_Boyd_061522/CAFI_stand_age/AK_PSP_increment_cores USU Mar 2013b.xlsx', sheet = 1, skip = 1) %>%
  filter(if_any(everything(), ~ !is.na(.)))

# fix the PSP
tree.2 <- trees %>%
  mutate(sub1 = substr(PSP, 1,3),
    subp = as.numeric(str_extract(sub1, "^\\d+")))

unique.subs <- subs %>%
  select(PLOT, SUBP)%>%
  distinct() %>%
  mutate(subp = row_number())

mel.age <- mel_age %>%
  rename(PLOT = site) %>%
  mutate(Spp = case_when(
    Sp_age_measured == "birch"  ~ "benteo",
    Sp_age_measured == "white spruce"  ~ "picgla",
    Sp_age_measured == "aspen"  ~ "poptre",
    Sp_age_measured ==  "black spruce"  ~ "picmar",
    Sp_age_measured ==  "poplar"  ~ "popbal",
    TRUE ~ 'NA'
  )) %>%
  select(-Sp_age_measured)

tree.test <- tree.2 %>%
  left_join(unique.subs) 

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

