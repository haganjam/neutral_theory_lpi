
# Title: Neutral theory as a null hypothesis for the Living Planet Index

# Project: Examine the Living Planet Index data

# load libaries
library(tidyverse)
library(here)

# load the LPI data
lpi_dat_raw <- read_csv(file = here("data/LPI_LPR2016data_public.csv"), 
                        na = c("", "NULL"),
                        col_types = cols(`2015` = col_double()))

# check which variables are included
names(lpi_dat_raw)

# check for any parsing problems
problems(lpi_dat_raw)

# check the data structure
str(lpi_dat_raw)

# View the data
View(lpi_dat_raw)


# subset the data based on meeting

# subset out the mammals
lpi_dat_raw$Class %>%
  unique()

lpi_dat_raw <- 
  lpi_dat_raw %>%
  filter(Class == "Mammalia")

# subset out populations with specific locations
lpi_dat_raw$Specific_location %>%
  unique()

lpi_dat_raw <- 
  lpi_dat_raw %>%
  filter(Specific_location == 1)

# subset out populations with actual count data
lpi_dat_raw$Units %>%
  unique()

lpi_dat_raw$Method %>%
  unique()

select(lpi_dat_raw, Units, Method) %>%
  View()

lpi_dat_raw %>%
  filter(grepl(pattern = "Indiv", x = Units) | grepl(pattern = "popula", x = Units)) %>%
  pull(Units) %>%
  unique()

# how do we choose which populations to include?




















                        