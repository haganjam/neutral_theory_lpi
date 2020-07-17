
# Title: Neutral theory as a null hypothesis for the Living Planet Index

# Project: Santini (2018): Global drivers of population abundance in terrestrial vertebrates

# test how to predict vertebrate population density from their model

sp_n <- 100

# population density estimates
dens <- rnorm(n = sp_n, mean = 150, sd = 40)
dens

# family
fam <- sample(c("bird", "mammal", "amphibian", "reptile"),
              size = sp_n, replace = TRUE)
fam

# species
spp <- sample(as.character(c(1:(sp_n/3) )), size = sp_n, replace = TRUE)
spp

dat <- 
  data.frame(id = 1:length(fam),
             family = fam)

dat

n_fam <- length(unique(dat$family))

spp <- 
  list("1" = as.character(1:10),
     "2" = as.character(20:30),
     "3" = as.character(40:50),
     "4" = as.character(70:90))

for (i in seq_along(1:n_fam)) {
  
  dat[dat$family == unique(dat$family)[i],]$species <- 
    sample(spp[[i]], size = nrow(dat[dat$family == unique(dat$family)[i],]), replace = TRUE)
  
}




# numeric predictor e.g. body size
bm <- 500 + ((-1)*dens) + rnorm(n = sp_n, mean = 0, sd = 30)
bm

# set up a location
loc <- sample(LETTERS, size = sp_n, replace = TRUE)
loc

# bind these into a dataframe
dat <- 
  data.frame(location = loc,
             family = fam,
             species = spp,
             body_mass = bm,
             density = dens)

length(unique(dat$species))

# load the lme4 library
library(lme4)
library(piecewiseSEM)

# fit a model to these data

lmm_1 <- lmer(formula = density ~ body_mass + (1|family/species) + (1|location),
              data = dat, REML = TRUE)

rsquared(lmm_1)

summary(lmm_1)

predict(lmm_1, type = "response")

ranef(lmm_1)





