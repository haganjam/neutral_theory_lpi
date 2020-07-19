
# Title: Neutral theory as a null hypothesis for the Living Planet Index

# Project: Santini (2018): Global drivers of population abundance in terrestrial vertebrates

# test how to predict vertebrate population density from their model

n_fam <- 4

n_spp <- 50

sim_dat <- data.frame(family = rep(LETTERS[1:n_fam], each = n_spp*n_fam),
           species = rep(as.character(1:50), times = n_fam))

# make the species' names unique
sim_dat$species <- paste(sim_dat$family, sim_dat$species, sep = "_")
sim_dat

# set up a vector of mean body mass values for each species
bm <- c(10, 50, 100, 200)

bm_out <- vector("list", length = length(bm))
for (i in 1:length(bm)) {
  
  bm_out[[i]] <- 
    rnorm(n = n_spp, mean = bm[i], sd = bm[i]/2) + rnorm(n = n_spp, mean = 0, sd = 5)
  
}

# add these body mass data to the sim_dat
sim_dat$body_mass <- unlist(bm_out)

sim_dat

# numeric predictor e.g. body size
sim_dat$density <- 1000 + ((-1)*sim_dat$body_mass) + rnorm(n = nrow(sim_dat), mean = 0, sd = 50)
sim_dat

# set up a location
n_loc <- 200

sim_dat$location <- sample(as.character(c(1:n_loc)), size = nrow(sim_dat), replace = TRUE)
sim_dat


# load the lme4 library
library(lme4)
library(piecewiseSEM)

# fit a model to these data

lmm_1 <- lmer(formula = density ~ body_mass + (1|family/species) + (1|location),
              data = sim_dat, REML = TRUE)

rsquared(lmm_1)

summary(lmm_1)

predict(lmm_1, type = "response")

ranef(lmm_1)

# given this model, how do we predict the values for a new point?

# location = 10

sim_dat[sim_dat$family == "A" & sim_dat$location == 179, ]

predict(lmm_1, data.frame(family = "A",
           location = "10",
           species = NA,
           body_mass = 20),
        type = "response")

# we probably want to run the model again but instead of having species in the model,
# just have the family in the model...

# seems reasonable.


