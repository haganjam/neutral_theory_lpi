
# Title: Neutral theory as a null hypothesis for the Living Planet Index

# Project: Random Poisson process and the Living Planet Index (Falko)

set.seed(41)
par(mai=c(0.8,0.8,0.2,0.2))

# This is the duration of the simulation, in year, to coincide with the LPI
years <- 1970:2020

# This is the lambda parameter in the Poisson distribution
Lambda <- 2

# This calls a blank plot, onto which each iterated simulation will be plotted
plot(0,0,type="n", las=1,xlab="Years", ylab ="Population",ylim=c(0,200), xlim=c(1970,2020))

N.start <- 100 # The starting population of the focal species

#Iterate the model 100 times for imaginary species
for (k in 1:100) {
  N <- rep(NA,length(years)); N[1] <- N.start
  for (i in 2:length(years)) {
    if(N[i-1] > 0) {
      
      # Population in time Nt+1 is equal to population Nt, minus deaths, plus births
      N[i] <- N[i-1] + (sample(c(-1,1),1)*rpois(1,Lambda))
    } else {
      N[i] <- 0 ; N[i-1] <- 0
    }
  }
  # Add lines to the plot
  lines(years,N,col=rgb(1,0,0,0.2))
}


# Starting population size at t=1
Nt1 <- 75

#Fluctuations
Fluct <- -50:50

# Population at time t=2
Nt2 <- Nt1 + Fluct

# Log response-ration (Lambda)
lam.75 <- log10(Nt2/Nt1)

# Make the plot
plot(Fluct,lam.75,las=1,xlab="Size of fluctuation (delta N)",
     ylab="Log response-ratio (lambda)", col=rgb(1,0,0,0.7))

# Repeat for starting population N=100
Nt1 <- 100
Nt2 <- Nt1 + Fluct
lam.100 <- log10(Nt2/Nt1)

points(Fluct,lam.100, col=rgb(0,0,1,0.7))


# Repeat for starting population N=150
Nt1 <- 150
Nt2 <- Nt1 + Fluct
lam.150 <- log10(Nt2/Nt1)
points(Fluct,lam.150, col=rgb(0,1,0,0.7))
# Add a legend
legend("topleft", pch=1, cex=0.7,
       col=c("green", "blue", "red"),
       c("Nt1 = 150", "Nt1 = 100", "Nt1 = 75"))



# see how this problem affects the LPI

set.seed(41)
years <- 1970:2020
N.start <- 100
Lambda <- 5
# Create a data matrix for the stochastic populations (adding two columns)
Pop.data <- matrix(NA, nrow=100, ncol=length(years)+2)
# The LPI requires a specific format. The first column is a ID vector of time-series
Pop.data[,1] <- 1:100

# The second column is the 'species names'.
#Here, assume that there are 100 unique species
Pop.data[,2] <- as.factor(1:100)
# Add columns names (note: LPI requires that years be preceded by an 'X")
colnames(Pop.data) <- (c("ID","Binomial",paste0("X",as.factor(years))))

# This is esential the same simulation as above
for (k in 1:100) {
  N <- rep(NA,length(years)); N[1] <- N.start
  for (i in 2:length(years)) {
    if(N[i-1] > 0) {
      # Population in time Nt+1 is equal to population Nt, minus deaths, plus births
      N[i] <- N[i-1] + (sample(c(-1,1),1)*rpois(1,Lambda))
    } else {
      N[i] <- 0 ; N[i-1] <- 0
    }
  }
  Pop.data[k,3:53] <- N
}

Pop.data

# Load the 'rlpi' package
library(rlpi)

# This is just an index vector for wich time-serie hould be included in calculations
# Here, we include all time series

# (Note: this command in needed when you want to use a subset of the data)
index_vector <- rep(TRUE, nrow(Pop.data))

#This creates an 'infile' for calcualting the LPI

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
     paste("N(start) = ", N.start, ", Poisson(lambda) = ", Lambda))




