
# Title: Neutral theory as a null hypothesis for the Living Planet Index

# Project: Examine the Living Planet Index data

# load libaries
library(tidyverse)
library(here)
library(vegan)

# load the LPI data
lpi_dat_raw <- read_csv(file = here("data/LPI_LPR2016data_public.csv"), 
                        na = c("", "NULL"),
                        col_types = cols(`2015` = col_double()))

lpi_dat_raw

# check which variables are included
names(lpi_dat_raw)

# check for any parsing problems
problems(lpi_dat_raw)

# check the data structure
str(lpi_dat_raw)

# View the data
View(lpi_dat_raw)

# remove years before 1970 and after 2012 because there is insufficient data (McCrae et al. 2017)
lpi_dat_raw <- 
  lpi_dat_raw %>%
  select(!starts_with(c("195", "196", "2013", "2014", "2015")))

# GAMs are used to extroplate time-series with > 6 data points between start and end years
# i.e. they don't extrapolate beyond range of the data

# Chain method is used for time-series with < 6 data points
# log-linear imputation when year(-1) = X, year(0) = unknown, year(+1) = known
# the unknown year is then imputed

# this makes it complicated to determine how this is calculated
# specifically, which series are used in the calculation for each year


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


# how many unique population time-series are there?
n_series <- 
  lpi_pop$ID %>%
  unique(.) %>%
  length(.)

# how many unique population time-series have more than 6 data points?
n_series_6 <- 
  lpi_pop %>%
  filter(!is.na(population_size)) %>%
  group_by(ID) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 6) %>%
  nrow()

n_series_6/n_series # <50% have more than 6 data points

# how many starting years are after 1970?
pop_id <- 
  lpi_pop %>%
  filter(!is.na(population_size)) %>%
  group_by(ID) %>%
  summarise(min_year = min(year),
            max_year = max(year),
            n = n(),
            .groups = "drop")


# gam method (definitely)
gam_id <- 
  pop_id %>%
  filter(n > 6) %>%
  pull(ID)


# create a variable with years between start and end date of each population
# 'Population time series with fewer than six data points or that resulted in poor GAM fit were modelled using the chain method [9]'
# from McCrae et al. (2017)

n0_nt <- 
  pop_id %>%
  select(min_year, max_year) %>%
  split(., pop_id$ID)

n0_nt <- 
  lapply(n0_nt, function(x) {
  tibble(year = seq(from = x$min_year, to = x$max_year, by = 1)) } ) %>%
  bind_rows(., .id = "ID") %>%
  mutate(ID = as.numeric(ID),
         year = as.character(year))

# join these data to the lpi pop data
n0_nt_dat <- 
  inner_join(n0_nt, lpi_pop, by = c("ID", "year"))

# in the LPI calculation, you need two consecutive years
# remove the first year for each of these time-series

n0_nt_dat <- 
  n0_nt_dat %>%
  group_by(ID) %>%
  filter(year != min(year)) %>%
  ungroup()

n0_nt_dat

# create a species x year matrix (columns = species, rows = years)
spp_year <- 
  split(n0_nt_dat, n0_nt_dat$year) %>%
  lapply(., function(x) {
    
    unique(x$Class)
    
  })

# create a vector of all species in the database
spp_list <- unique(n0_nt_dat$Class)

spp_year <- lapply(spp_year, function(x) { 
  
  tibble(species = spp_list,
         pa = if_else(spp_list %in% x, 1, 0))
  
  } ) %>%
  bind_rows(., .id = "year")

spp_year <- 
  spp_year %>%
  pivot_wider(names_from = species,
              values_from = pa)

# run an nmds on these data to see how the composition has changed through time
spp_year

nmds.1 <- 
  metaMDS(select(spp_year, -year), 
          distance = "jaccard", k = 2, try = 20, trymax = 20, 
          engine = c("monoMDS"), autotransform = FALSE, 
          wascores = TRUE, expand = TRUE, 
          trace = 1, plot = FALSE)

stressplot(nmds.1)

nmds.1_raw <- as_tibble(nmds.1$points)

# add the years to this dataframe
nmds.1_raw <- 
  nmds.1_raw %>%
  mutate(year = as.numeric(spp_year$year))

ggplot(data = nmds.1_raw,
       mapping = aes(x = year, y = MDS1, colour = year)) +
  geom_point() +
  scale_colour_viridis_c() +
  theme_classic()

# use multidimensional scaling (MDS) i.e. principle coordinates analysis
mds.1 <- cmdscale(vegdist(select(spp_year, -year), method = "jaccard"), 
                  k = 2, eig = T, add = T )

# calculate the variance explained by the mds.1
round(mds.1$eig*100/sum(mds.1$eig), 1)

# check the correlation between observed differences and actual dissimilarities
round(cor(vegdist(select(spp_year, -year), method = "jaccard"), 
          dist(mds.1$points), method = "spearman"), 3)

# create a dataframe for plotting
mds.1_raw <- 
  tibble(mds_1 = mds.1$points[, 1],
         mds_2 = mds.1$points[, 2]) %>%
  mutate(year = as.numeric(spp_year$year))

ggplot(data = mds.1_raw,
       mapping = aes(x = year, y = mds_1, colour = year)) +
  geom_jitter(size = 2.5, alpha = 0.9) +
  scale_colour_viridis_c() +
  ylab("MDS1 (43 %)") +
  theme_classic()


# analyse the number of time-series in different classes

prop_class <- 
  n0_nt_dat %>%
  group_by(year, Class) %>%
  summarise(n_class = n(), .groups = "drop") %>%
  group_by(year) %>%
  mutate(n_total = sum(n_class)) %>%
  ungroup() %>%
  mutate(class_proportion = n_class/n_total) %>%
  mutate(year = as.numeric(year))

ggplot(data = prop_class,
       mapping = aes(x = year, y = class_proportion, colour = Class)) +
  geom_point() +
  geom_smooth(method = "lm", size = 0.1, se = TRUE, alpha = 0.1) +
  scale_colour_viridis_d() +
  theme_classic()


# analyse the number of time-series in different systems

prop_syst <- 
  n0_nt_dat %>%
  group_by(year, System) %>%
  summarise(n_syst = n(), .groups = "drop") %>%
  group_by(year) %>%
  mutate(n_total = sum(n_syst)) %>%
  ungroup() %>%
  mutate(syst_proportion = n_syst/n_total) %>%
  mutate(year = as.numeric(year))

ggplot(data = prop_syst,
       mapping = aes(x = year, y = syst_proportion, colour = System)) +
  geom_point() +
  geom_smooth(method = "lm", size = 0.1, se = TRUE, alpha = 0.1) +
  scale_colour_viridis_d() +
  theme_classic()


# extract the starting population size for each species
start_pop <- 
  lpi_pop %>%
  filter(!is.na(population_size)) %>%
  group_by(ID) %>%
  filter(year == min(year)) %>%
  ungroup()

# check the distribution of starting population sizes
summary(start_pop$population_size)

# summarise the starting population data
start_pop <- 
  start_pop %>%
  mutate(population_size = log10(1 + population_size)) %>%
  group_by(Class, year) %>%
  summarise(population_size_m = mean(population_size, na.rm = TRUE),
            population_size_sd = sd(population_size, na.rm = TRUE),
            n = n(),
            .groups = "drop") %>%
  mutate(population_size_se = (population_size_sd/sqrt(n)),
         year = as.numeric(year)) 

# how have starting population sizes changed through time?
ggplot(data = start_pop,
       mapping = aes(x = as.numeric(year), y = log10(1 + population_size_m))) +
  geom_jitter(mapping = aes(),
              alpha = 0.25, shape = 16, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_viridis_d() +
  xlab("year") +
  theme_classic()


# population sizes through time for the different groups
pop_group <- 
  lpi_pop %>%
  filter(!is.na(population_size)) %>%
  mutate(population_size = log10(1 + population_size)) %>%
  group_by(Class, year) %>%
  summarise(population_size_m = mean(population_size, na.rm = TRUE),
            population_size_sd = sd(population_size, na.rm = TRUE),
            n = n(),
            .groups = "drop") %>%
  mutate(population_size_se = (population_size_sd/sqrt(n)),
         year = as.numeric(year)) 
  
ggplot(data = pop_group,
       mapping = aes(x = year, y = population_size_m)) +
  geom_errorbar(mapping = aes(ymin = population_size_m-population_size_se,
                              ymax = population_size_m+population_size_se,
                              colour = Class),
                width = 0.1) +
  geom_point(mapping = aes(colour = Class),
             alpha = 1, shape = 16, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_viridis_d() +
  xlab("year") +
  theme_classic()


# population sizes through time for all groups combined
pop_all <- 
  lpi_pop %>%
  filter(!is.na(population_size)) %>%
  mutate(population_size = log10(1 + population_size)) %>%
  group_by(year) %>%
  summarise(population_size_m = mean(population_size, na.rm = TRUE),
            population_size_sd = sd(population_size, na.rm = TRUE),
            n = n(),
            .groups = "drop") %>%
  mutate(population_size_se = (population_size_sd/sqrt(n)),
         year = as.numeric(year)) 

pop_all

p_all <- 
  ggplot(data = pop_all,
       mapping = aes(x = year, y = population_size_m)) +
  geom_errorbar(mapping = aes(ymin = population_size_m-population_size_se,
                              ymax = population_size_m+population_size_se),
                width = 0.1) +
  geom_point(alpha = 1, shape = 16, size = 2) +
  geom_smooth(size = 0.1, colour = "black") +
  scale_colour_viridis_d() +
  xlab("year") +
  ylab("mean +- se population size") +
  theme_classic()

p_all


# how do we choose which populations to include?




















                        