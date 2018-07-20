#-------------------------------------------------------------------------------
# IBPNEW 2012
# Code to generate the .dat file for input into the ADMB model
# Uses Lowestoft format input files
# David Miller, some code adapted from FLR-project
# Date: 3 Oct 2012
#-------------------------------------------------------------------------------
# FLR
# install.packages(repos="http://flr-project.org/R")
## Source methods/functions
source(paste(codePath,"nsea_functions.r",sep=""))

### ------------------------------------------------------------------------------------------------------
###  1. Settings
### ------------------------------------------------------------------------------------------------------

## Assessment settings
## Stock name
Units               <- c("tonnes","thousands","kg")     # totals (product of the next two), numbers, weights-at-age
maxA                <- 10
pGrp                <- T
minFbar             <- 2
maxFbar             <- 6

### ------------------------------------------------------------------------------------------------------
###   2. Read and process assessment input data
### ------------------------------------------------------------------------------------------------------

## Read stock data
stock               <- readFLStock(paste(dataPath,"index.txt", sep=""))
units(stock)[1:17]  <- as.list(c(rep(Units,4), "NA", "NA", "f", "NA", "NA"))
stkWtPgrp           <- stock.wt(stock)[maxA,]
if(pGrp)
  stock             <- setPlusGroup(stock, plusgroup=maxA) else stock <- trim(stock, age=1:maxA)
stock.wt(stock)[maxA,] <- stkWtPgrp
range(stock)[c("minfbar","maxfbar")] <- c(minFbar,maxFbar)
# Number of years
years               <- as.numeric(range(stock)["minyear"]:range(stock)["maxyear"])
numYr               <- as.numeric(range(stock)["maxyear"]-range(stock)["minyear"]+1)
# Number of ages
numAges             <- length(1:maxA)
startyr             <- range(stock)[["minyear"]]## Read index data

#- Setup indices
indices             <- readFLIndices(paste(dataPath,"fleet.txt", sep=""), na.strings="-1")
indices[[1]]@type   <- "number"
indices[[2]]@type   <- "number"
indices[[3]]@type   <- "biomass"
indices[[4]]@type   <- "number"
indices[[5]]@type   <- "biomass"
for(i in 1:length(indices))
  indices[[i]]@name <- names(indices)[i]

discards.wt(stock) <- catch.wt(stock) <- landings.wt(stock)

# Correct for SOP in landings.n
landings.n(stock)@.Data[which(landings.n(stock)==-1)] <- NA
landings.n(stock)   <- sweep(landings.n(stock),2,computeLandings(stock)/landings(stock),"/")
landings.n(stock)@.Data[which(is.na(landings.n(stock)))] <- -1
catch.n(stock)      <- landings.n(stock)
