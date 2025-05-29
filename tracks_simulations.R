#### SIMULATION ####

# Setting working directory
setwd("/import/ecoc9/data-jeltsch/arceguillen/")

##### Packages needed #####

library(raster)
library(CircStats)
library(geoR)
library(tictoc)
library(foreach)
library(doParallel)
library(circular)
library(copula)
library(cylcop)
library(distr)
library(fmesher)
library(mgcv)
library(ggplot2)
library(sp)
library(sf)
library(fields)

# Parallel kernels specification
library(parallel)
cl <- makeCluster(80)

# Source function to simulate data
source("fun_complete_gam.R")

##### GENERATION LANDSCAPES X1, X2 and cen #####

# Gaussian random field with an exponential covariance function

# The marginal variance equal to 1

# Defining extend
xextent <- 100
yextent <- 100
grid_resolution <- 0.1

# Creating grid list for Gaussian random field
grid <- list(x=seq(-xextent, xextent, grid_resolution),
             y= seq(-yextent, yextent, grid_resolution))

# Generating Gaussian random field

# X1 landscape with spatial range equal to 25 space units #
obj1 <- circulantEmbeddingSetup(grid, Covariance="Exponential", aRange= 25)

# Generating landscape
set.seed(1789)
gaussian_random_field1 <- circulantEmbedding(obj1)

# create raster object
landscape_1 <- raster(gaussian_random_field1,
                     xmn=-xextent-0.5*grid_resolution,
                     xmx=xextent+0.5*grid_resolution,
                     ymn=-yextent-0.5*grid_resolution,
                     ymx=yextent+0.5*grid_resolution)

# X2 landscape with spatial range equal to 15 space units #
obj2 <- circulantEmbeddingSetup(grid, Covariance="Exponential", aRange= 15)

# Generating landscape
set.seed(189)
gaussian_random_field2 <- circulantEmbedding(obj2)

# create raster object
landscape_2 <- raster(gaussian_random_field2,
                     xmn=-xextent-0.5*grid_resolution,
                     xmx=xextent+0.5*grid_resolution,
                     ymn=-yextent-0.5*grid_resolution,
                     ymx=yextent+0.5*grid_resolution)

names(landscape_1) = "x1"
names(landscape_2) = "x2"
cor(values(landscape_1), values(landscape_2))


# Creating centralizing tendency cen
# Defining starting point
cen <- c(0, 0)
cenmap <- distanceFromPoints(landscape_1, cen)
names(cenmap) <- "cen"

#### SCENARIOS ####
scenarios <- c("gamma_mises", "copula", "bimodal", "weibull", "copula_2")

#### GENERATING ANIMAL TRACKS FOR DIFFERENT SCENARIOS ####
for (case in scenarios) {
  # time series length
  Ts <- 1020
  
  # Making arbitrary times
  time.begin <- as.POSIXct("2017-11-25 12:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Berlin")
  dates <- time.begin + seq(0, (Ts - 1) * 360, 360)
  
  parameter <- list()
  
  # Defining true spatial effects values
  parameter$beta <- c(1.5, -1, -0.02)
  
  ##### SIMULATION #####
  
  # Joining landscapes to a single raster object
  landscape <- stack(landscape_1, landscape_2, cenmap)
  
  # Export objects needed for the parallel cluster
  clusterExport(cl = cl, varlist=c("scenarios", "landscape_1", "cen", "case", "landscape_2", "landscape", "turnAngle", "simdata_generic", "Ts", "time.begin", "parameter", "dates"),  envir = .GlobalEnv)
  
  iteration_list = as.list(1:100)
  
  # Creating tracks
  my_data <- parLapply(cl, iteration_list, function(i) {
    library(raster)
    library(CircStats)
    library(geoR)
    library(tictoc)
    library(foreach)
    library(doParallel)
    library(circular)
    library(copula)
    library(cylcop)
    library(distr)
    library(parallel)
    
    result <- simdata_generic(
      parameter = parameter, case = case, center = cen, Ts = Ts, seed = i,
      vzero = 1e-6, landscape = landscape
    )
    result <- result[21:Ts, ] # first 20 steps for initialisation
    result$t <- dates[21:Ts]
    return(result)
  })
  
  setwd("/import/ecoc9/data-jeltsch/arceguillen/")
  
  save(list = c("my_data", "landscape"), file = paste("gam_simulation_", case, ".RData", sep = ""))
}

