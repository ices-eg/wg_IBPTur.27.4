#-------------------------------------------------------------------------------
# IBPNEW 2012
# Code to generate the .dat file for input into the ADMB model
# Uses Lowestoft format input files
# David Miller, some code adapted from FLR-project
# Date: 3 Oct 2012
#-------------------------------------------------------------------------------
rm(list=ls())

# FLR
# install.packages(repos="http://flr-project.org/R")
library(FLCore);library(mgcv)
library(FLAssess);
library(stockassessment)
library(FLEDA); library(splines); 
library(scales); library(gplots);library(grid); library(gridExtra); library(latticeExtra)
library(sas7bdat)
library(TMB); library(FLSAM)

# Set paths to folders
Path      <- "D:/Repository/Turbot/assessment runs/"
dataPath  <- paste(Path,"Lowestoft files/",sep="")
outPath   <- paste(Path,"trial_runs_2017/Output/",sep="")
codePath  <- paste(Path,"Trial_runs_2017/",sep="")

## Source methods/functions
source(paste(codePath,"03a_setupStockIndices.r",sep=""))


run       <- "addBTSDG"
sens      <- "nopg_BTS_DG"

### ------------------------------------------------------------------------------------------------------
###   2. Read and process assessment input data
### ------------------------------------------------------------------------------------------------------

indices             <- FLIndices(list(window(trim(indices[[1]],age=1:6),start=2004),indices[[3]],indices[[16]]))
indices[["BTS_DG"]]@type <- "number"
indices[["BTS_DG"]] <- trim(indices[["BTS_DG"]],age=1:4)
### ------------------------------------------------------------------------------------------------------
###   3. Setup data structure for SAM assessment
### ------------------------------------------------------------------------------------------------------

TUR                 <- stock
TUR.tun             <- indices
TUR.ctrl            <- FLSAM.control(TUR,TUR.tun)

TUR.ctrl@states["catch",]                   <- c(0:6,rep(7,3))
TUR.ctrl@cor.F                              <- 2
TUR.ctrl@catchabilities["SNS",ac(1:6)]      <- c(0:2,rep(3,3))          + 101
TUR.ctrl@catchabilities["BTS_DG",ac(1:4)]   <- c(0,0,1,2)               + 201
TUR.ctrl@catchabilities["NL_LPUE",ac(1)]    <- 0                        + 301
TUR.ctrl@f.vars["catch",]                   <- c(0,1,2,2,3,3,3,4,4,4)
TUR.ctrl@logN.vars[]                        <- c(0,rep(1,9))
TUR.ctrl@obs.vars["catch",]                 <- c(0,0,1,1,1,2,2,2,3,3)   + 101
TUR.ctrl@obs.vars["SNS",ac(1:6)]            <- c(0,0,1,1,2,2)           + 201
TUR.ctrl@obs.vars["BTS_DG",ac(1:4)]         <- c(0,0,1,1)               + 301
TUR.ctrl@obs.vars["NL_LPUE",ac(1)]          <- 0                        + 401
TUR.ctrl@cor.obs[]                          <- NA
TUR.ctrl@cor.obs["SNS",1:5]                 <- c(0,0,rep(1,3))          + 101
TUR.ctrl@cor.obs["BTS_DG",1:3]              <- rep(0,3)                 + 201
TUR.ctrl@cor.obs.Flag[c(2,4)]               <- af("AR")
TUR.ctrl@biomassTreat[3]                    <- 2
TUR.ctrl                                    <- update(TUR.ctrl)
### ------------------------------------------------------------------------------------------------------
###   4. Run assessment
### ------------------------------------------------------------------------------------------------------

TUR.sam             <- FLSAM(TUR,TUR.tun,TUR.ctrl)
TUR.ctrl@residuals  <- FALSE; TUR.sam@control@residuals <- FALSE
TUR.retro           <- retro(TUR,TUR.tun,TUR.ctrl,retro=7,base.assess=TUR.sam)

### ------------------------------------------------------------------------------------------------------
###   5. Diagnostics
### ------------------------------------------------------------------------------------------------------
source(file.path(codePath,"03b_runDiagnostics.r"))

