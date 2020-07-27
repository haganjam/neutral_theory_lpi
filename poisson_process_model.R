
# Title: Neutral theory as a null hypothesis for the Living Planet Index

# Project: Random Poisson process and the Living Planet Index (Falko)


# how does the LPI respond to different randomly fluctuating populations

# Poisson process simulation

set.seed(41)
years <- 1970:2020
N.start <- 100
Lambda <- 5

# Create a data matrix for the stochastic populations (adding two columns)
Pois.data <- matrix(NA, nrow=100, ncol=length(years)+2)

# The LPI requires a specific format. The first column is a ID vector of time-series
Pois.data[,1] <- 1:100

# The second column is the 'species names'.

# Here, assume that there are 100 unique species
Pois.data[,2] <- as.factor(1:100)

# Add columns names (note: LPI requires that years be preceded by an 'X")
colnames(Pois.data) <- (c("ID","Binomial",paste0("X",as.factor(years))))

# This is esential the same simulation as above
for (k in 1:N.start) {
  N <- rep(NA,length(years)); N[1] <- N.start
  for (i in 2:length(years)) {
    if(N[i-1] > 0) {
      # Population in time Nt+1 is equal to population Nt, minus deaths, plus births
      N[i] <- N[i-1] + (sample(c(-1,1),1)*rpois(1,Lambda))
    } else {
      N[i] <- 0 ; N[i-1] <- 0
    }
  }
  Pois.data[k,3:53] <- N
}

Pois.data


# Normal distribution simulation

set.seed(41)
years <- 1970:2020
N.start <- 100
mean_fluc <- 0

# Create a data matrix for the stochastic populations (adding two columns)
norm.data <- matrix(NA, nrow=100, ncol=length(years)+2)

# The LPI requires a specific format. The first column is a ID vector of time-series
norm.data[,1] <- 1:100

# The second column is the 'species names'.

# Here, assume that there are 100 unique species
norm.data[,2] <- as.factor(1:100)

# Add columns names (note: LPI requires that years be preceded by an 'X")
colnames(norm.data) <- (c("ID","Binomial",paste0("X",as.factor(years))))

# This is esential the same simulation as above
for (k in 1:nrow(norm.data)) {
  N <- rep(NA,length(years)); N[1] <- N.start
  for (i in 2:length(years)) {
    if(N[i-1] > 0) {
      # Population in time Nt+1 is equal to population Nt, minus deaths, plus births
      N[i] <- round((N[i-1] + rnorm(n = 1, mean = mean_fluc, sd = 5)), 0)
    } else {
      N[i] <- 0 ; N[i-1] <- 0
    }
  }
  norm.data[k,3:53] <- N
}

norm.data

# Cyclic population simulation

# define generate_species function from Fournier et al. (2019)

# mu = mean population
# rho = autocorrelation 
# CV = Coefficient of Variation

# spp = number of species

generate_species <- function(nyears = 100, mu = 1000, rho = 0.5, CV. = 0.2, spp. = 20) {
  
  sim.years <- nyears+10
  
  ln_SD <- sqrt(log(CV.^2+1))  #SD parameter of lognormal distribution with desired variability
  
  ln_SD_e <- sqrt(ln_SD^2/(1-rho^2))  #SD for first value
  
  time_series_temp <- matrix(rnorm(mean = 0, sd = 1, n = sim.years*spp.),
                             ncol = spp., nrow = sim.years)
  
  time_series <- matrix(nrow = sim.years, ncol = spp.)
  
  time_series[1,] <- ln_SD_e*time_series_temp[1,]
  
  for (i in 2:(sim.years)) {
    time_series[i,] <- rho*time_series[(i-1),]+time_series_temp[i,]*ln_SD
  }
  
  time_series <- time_series + log(mu) - 0.5*ln_SD_e^2
  
  time_series <- exp(time_series)
  
  return(time_series[-c(1:10),])
}

# Create a data matrix for the stochastic populations (adding two columns)
cyclic.data <- matrix(NA, nrow=100, ncol=length(years)+2)

# The LPI requires a specific format. The first column is a ID vector of time-series
cyclic.data[,1] <- 1:100

# The second column is the 'species names'.

# Here, assume that there are 100 unique species
cyclic.data[,2] <- as.factor(1:100)

# Add columns names (note: LPI requires that years be preceded by an 'X")
colnames(cyclic.data) <- (c("ID","Binomial",paste0("X",as.factor(years))))


sim_out <- t(generate_species(nyears = length(1970:2020), 
            mu = 100, rho = 0.5, CV. = 0.2, spp. = 100) )

for (k in 1:100) {
  
  cyclic.data[k , 3:53] <- sim_out[k, ]
  
}

cyclic.data


# Load the 'rlpi' package
library(rlpi)

# select a data set to use
Pop.data <- cyclic.data

# This is just an index vector for wich time-series should be included in calculations
# Here, we include all time series

# (Note: this command in needed when you want to use a subset of the data)
index_vector <- rep(TRUE, nrow(Pop.data))

#This creates an 'infile' for calcualting the LPI

# Note: you need a folder names 'LPI_files' in your working directory
library(here)

if(! dir.exists(here("LPI_files"))){
  dir.create(here("LPI_files"))
}

sim_infile <- create_infile(as.data.frame(Pop.data), start_col_name="X1970",
                            end_col_name="X2020", index_vector=index_vector, name="LPI_files/lpi_mam")

# Calclualte LPI using this infile, for the period 1970 to 2019 with 100 bootstraps.
sim_lpi <- LPIMain(sim_infile, REF_YEAR = 1970, PLOT_MAX = 2019,
                   BOOT_STRAP_SIZE = 10, VERBOSE=FALSE, plot_lpi=FALSE)


# Plot the output
par(mai=c(0.8,0.8,0.2,0.2))
plot(0,0,type="n",las=1,xlim=c(1970,2020),
     ylim=c(0.8,1.2),ylab="LPI (1970 = 1)",xlab="Year")

# Add baseline
abline(h=1,col="grey")

# Include error bars as polygons
polygon(c(seq(1970,2020),seq(2020,1970)),
        c(sim_lpi$CI_low,rev(sim_lpi$CI_high)),col=rgb(0,0.5,0,0.75),border=NA)

# Add mean line
lines(c(1970:2020),sim_lpi$LPI_final,col="white",lwd=2)

# Include text with simulation parameters
text(1995, 1.15,
     paste("N(start) = ", N.start, ", Poisson(lambda) = ", Lambda))




