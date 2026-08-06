## Appending the 2026 Data!

# Read in the 2026 data
library(readxl)
library(dplyr)
library(purrr)
library(stringr)
library(readr)
library(tidyr)


# 1. Get all matching files (ignore.case handles All/all/ALL)
files <- list.files("/Users/olhajek/Desktop/RSN/RSN_proj/Data/2026 Data/Completed Sites",
                    pattern = "^2026_.*_Trees_All\\.xlsx$",
                    full.names = TRUE,
                    ignore.case = TRUE)

# 2. Pull out the "XXX" portion, case-insensitively
site_names <- str_extract(basename(files),
                          regex("(?<=^2026_)[^_]+(?=_Trees_All\\.xlsx$)", ignore_case = TRUE))

# 3. Read every file with your function, naming each result by its site code

# create a function to read in the data
read_tree_sheet <- function(path) {
  
  # Read the ENTIRE sheet as raw text, no column names yet
  # (col_types = "text" also avoids the "Expecting numeric... got a date" warnings,
  # since nothing is being type-guessed at this stage)
  raw <- read_excel(path, col_names = FALSE, col_types = "text")
  
  row1 <- trimws(as.character(raw[1, ]))
  row2 <- trimws(as.character(raw[2, ]))
  
  col_names <- ifelse(is.na(row2) | row2 == "",
                      row1,
                      paste0(row1, "_", row2))
  
  # Make any duplicate column names unique (e.g., "NOTES", "NOTES" -> "NOTES", "NOTES_1")
  col_names <- make.unique(col_names, sep = "_")
  
  # Drop the two header rows; the rest is data
  data <- raw[-c(1, 2), ]
  names(data) <- col_names
  
  
  # Everything was read as text above -- convert numeric-looking columns back
  data <- type_convert(data, guess_integer = FALSE)
  
  data
}

# drop FP5D because format is different
# Drop FP5D (case-insensitive, in case it's sometimes "fp5d" or "Fp5d")
keep <- toupper(site_names) != "FP5D"
files <- files[keep]
site_names <- site_names[keep]

site_list <- files %>%
  set_names(site_names) %>%
  map(read_tree_sheet)


# Create a dataframe iwht a column lookup for this, will use this to merge 
column_name_lookup <- imap_dfr(site_list, ~ tibble(site = .y, raw_name = names(.x))) %>%
  mutate(clean_name = NA_character_)

# Read in FP5D and merge that here
FP5D <- read_excel('/Users/olhajek/Desktop/RSN/RSN_proj/Data/2026 Data/Completed Sites/2026_FP5D_Trees_all.xlsx') %>%
  # I think will want to coalesce the notes for FP5D
  mutate(NOTES = coalesce(...15, ...16,...17)) %>%
  rename(DBH = "2026.0")
  
glimpse(FP5D)
fp5d_colnames <- tibble(
  site = "FP5D",
  raw_name = colnames(FP5D),
  clean_name = NA_character_)

# Join FP5D into the column lookup
column_lookup <-rbind(column_name_lookup, fp5d_colnames)


# Export the column lookup and will manually add it 
#write.csv(column_lookup, '/Users/olhajek/Desktop/RSN/RSN_proj/Data/2026 Data/Data_checks/column_lookup_2026_data.csv', row.names = F)

# read in the CSV
match <- read.csv('/Users/olhajek/Desktop/RSN/RSN_proj/Data/2026 Data/Data_checks/column_lookup_2026_data.csv') %>%
  filter(!is.na(clean_name) & trimws(clean_name) != "")

# fix that NA name in UP4d
names(site_list$UP4D)
names(site_list$UP4D)[is.na(names(site_list$UP4D))] <- "NOTES"

names(site_list$BDm1)

site_list_clean <- imap(site_list, function(df, site_id) {
  
  site_lookup <- match %>%
    filter(site == site_id, raw_name %in% names(df))
  
  df %>%
    select(all_of(site_lookup$raw_name)) %>%
    rename(!!!setNames(as.list(site_lookup$raw_name), site_lookup$clean_name))
})

# join the site list
# lots of inconsistencies and errors in the data
map(site_list_clean, ~sapply(.x, class))

# going to convert everything to character then individually convert back
# combining here
all_sites <- map(site_list_clean, ~ mutate(.x, across(everything(), as.character))) %>%
  bind_rows(.id = "site")

# going to add in FP5D now!

site_lookup <- match %>%
  filter(site == "FP5D")

fp5d.2 <- FP5D %>%
  select(all_of(site_lookup$raw_name)) %>%
  rename(!!!setNames(as.list(site_lookup$raw_name), site_lookup$clean_name)) %>%
  mutate(across(everything(), as.character))
str(fp5d.2)

all_site2 <- all_sites %>%
  bind_rows(fp5d.2)

## Make sure to remove any totally blank rows, site, add in the same rows that are present in the raw data set
glimpse(all_site2)

unique(all_site2$SITE)
unique(all_site2$PLOT)
unique(all_site2$SPECIES)

all <- all_site2 %>%
  # remove site column 
  select(-site)%>%
  #remove any all blank rows
  filter(!if_all(everything(), is.na)) %>% 
  # remove these two rows because they are duplicates and will be corrected in the PLOT line
  filter(is.na(DBH) | DBH != "Entered in plot 2") %>%
  filter(is.na(DBH) | DBH != "this is now tree 357") %>%
  # start to clean up the data frame
  mutate(SITE = case_when(
    SITE == "BDMI" ~ "BDM1",
    SITE == "GAM4" ~ "GSM4",
    SITE == "Up4b" ~ "UP4B",
    TRUE ~ SITE
  ), 
  PLOT = case_when(
    is.na(PLOT) & SITE == "GSM2" ~ "12", 
    is.na(PLOT) & SITE == "UP4B" ~ "12",
    SITE == "UP4C" & TREE == "30" & PLOT == "2" ~ "3",
    SITE == "UP4C" & TREE == "31" & PLOT == "2" ~ "3",
    TRUE ~PLOT
  ), 
  TREE = ifelse(is.na(TREE), "2579", TREE), 
  # NOTE THERE IS an 802a and 801A in UP4B. neither of these show up before, but just adding them here now
  TREE = as.integer(str_remove(TREE, "[Aa]$")), 
  SPECIES = case_when(
    SPECIES == "Picmar" ~"PICMAR",
    SPECIES == "PICMA" ~"PICMAR",
    # assuming these trees are PICMAR; slight chance could be LARLAR, but likely PICMAR and these are NEW
    SPECIES == "UP4C" ~ "PICMAR", 
    is.na(SPECIES) ~ "PICMAR",
    TRUE ~ SPECIES), 
  DBH = case_when(
    str_starts(DBH, "-7") ~ "-7777",
    str_starts(DBH, "-8") ~ "-8888",
    str_starts(DBH, "-9") ~ "-9999",
    TRUE ~ DBH))

  mutate(UniqueID = paste(SITE,PLOT,TREE, sep="_"), 
         Date = NA_character_, DA) %>%
    # Drop these two rows because they were entered wrongly and they've been updated in the previous section
    filter(UniqueID %!in% c("UP4C_2_30", "UP4C_2_31"))

unique(all$SITE)
unique(all$PLOT)
unique(all$SPECIES)


# NOTES
# TREE 357 is the old 26 in UP4A Plot 1 - this doesn't seem to be updated in the final script either



