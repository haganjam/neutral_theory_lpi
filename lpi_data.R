
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


# determine how the population data has changed through time
# make a year and population size column

lpi_dat_raw %>%
  select(starts_with(c("19", "20"))) %>%
  names()

names(lpi_dat_raw)

lpi_pop <- 
  lpi_dat_raw %>%
  pivot_longer(cols = starts_with(c("19", "20")),
               names_to = "year",
               values_to = "population_size")


# extract the starting population size for each species
start_pop <- 
  lpi_pop %>%
  filter(!is.na(population_size)) %>%
  group_by(ID) %>%
  filter(year == min(year)) %>%
  ungroup()

summary(start_pop$population_size)


ggplot(data = start_pop,
       mapping = aes(x = as.numeric(year), y = log10(1 + population_size))) +
  geom_jitter(mapping = aes(colour = Class),
              alpha = 0.25, shape = 16, size = 2) +
  geom_smooth(method = "lm", colour = "black", size = 0.5) +
  scale_colour_viridis_d() +
  theme_classic()



# how do we choose which populations to include?




















                        