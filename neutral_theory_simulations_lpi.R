
# Title: Neutral theory as a null hypothesis for the Living Planet Index

# Project: Neutral theory simulations (Falko)

set.seed(41)
par(mai=c(0.8,0.8,0.2,0.2))

# This is the duration of the simulation, in year, to coincide with the LPI
years <- 1970:2020

# This calls a blank plot, onto which each iterated simulation will be plotted
plot(0,0,type="n", las=1,xlab="Years", ylab ="Population",ylim=c(0,200), xlim=c(1970,2020))
J <- 500 # Total number of ecological equivalent individuals
lifespan <- 5 # Lifespan of the focal species
N.start <- 100 # The starting population of the focal species

# Iterate the model 100 times
for (k in 1:100) {
  N <- rep(NA,length(years)); N[1] <- N.start
  for (i in 2:length(years)) {
    # This is the number of birth-deaths 'events' per year in the focal population
    # It means that their is complete population turnover at the scale of lifespans
    # We round up fraction
    events <- ceiling(N[i-1] *(1/lifespan))
    # The death process is proportional to the realtive abundance of the focal population
    death <- sample(c(1,0),events,prob=c(N[i-1]/J, (1-N[i-1]/J)),replace=TRUE)
    # The birth process is proportional to the realtive abundance of the focal population
    birth <- sample(c(1,0),events,prob=c(N[i-1]/J, (1-N[i-1]/J)),replace=TRUE)
    # Population in time Nt+1 is equal to population Nt, minus deaths, plus births
    N[i] <- N[i-1] - sum(death) + sum(birth)
  }
  # Add lines to the plot
  lines(years,N,col=rgb(1,0,0,0.2))
}


set.seed(41)
years <- 1970:2020
N.start <- 100
J <- 500
lifespan <- 5

# Create a data matrix for the stochastic populations (adding two columns)
Pop.data <- matrix(NA, nrow=100, ncol=length(years)+2)

# The LPI requires a specific format. The first column is a ID vector of time-series
Pop.data[,1] <- 1:100

# The second column is the 'species names'.
# Here, assume that there are 100 unique species
Pop.data[,2] <- as.factor(1:100)

# Add columns names (note: LPI requires that years be preceded by an 'X")
colnames(Pop.data) <- (c("ID","Binomial",paste0("X",as.factor(years))))
# This is esential the same simulation as above
for (k in 1:100) {
  N <- rep(NA,length(years)); N[1] <- N.start
  for (i in 2:length(years)) {
    events <- ceiling(N[i-1] *(1/lifespan))
    death <- sample(c(1,0),events,prob=c(N[i-1]/J, (1-N[i-1]/J)),replace=TRUE)
    birth <- sample(c(1,0),events,prob=c(N[i-1]/J, (1-N[i-1]/J)),replace=TRUE)
    N[i] <- N[i-1] - sum(death) + sum(birth)
  }
  Pop.data[k,3:53] <- N
}

# Load the 'rlpi' package
library(rlpi)
# This is just an index vector for wich time-serie hould be included in calculations
# Here, we include all time series

# (Note: this command in needed when you want to use a subset of the data)
index_vector <- rep(TRUE, nrow(Pop.data))

# This creates an 'infile' for calcualting the LPI
# Note: you need a folder names 'LPI_files' in your working directory
sim_infile <- create_infile(as.data.frame(Pop.data), start_col_name="X1970",
                            end_col_name="X2020", index_vector=index_vector, name="LPI_files/lpi_mam")

# Calclualte LPI using this infile, for the period 1970 to 2019 with 100 bootstraps.
sim_lpi <- LPIMain(sim_infile, REF_YEAR = 1970, PLOT_MAX = 2019,
                   BOOT_STRAP_SIZE = 100, VERBOSE=FALSE, plot_lpi=FALSE)

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
     paste("N(start) = ", N.start, ", J = ", J, ",\nLifespan = ", lifespan))









    