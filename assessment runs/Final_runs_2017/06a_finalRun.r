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

run       <- "final"
sens      <- "final"

### ------------------------------------------------------------------------------------------------------
###   2. Read and process assessment input data
### ------------------------------------------------------------------------------------------------------

indices             <- FLIndices(list(window(trim(indices[[1]],age=1:6),start=2004),window(indices[[2]],start=1991),indices[[9]]))
indices[["Dutch_BT2_LPUE_ModelD"]]@type <- "biomass"
names(indices)      <- c("SNS","BTS-ISIS","NL_LPUE")

### ------------------------------------------------------------------------------------------------------
###   3. Setup data structure for SAM assessment
### ------------------------------------------------------------------------------------------------------

TUR                 <- window(stock,start=1981)
TUR@catch.n[,ac(2000:2002)] <- -1; TUR@landings.n[] <- TUR@catch.n
TUR                 <- setPlusGroup(TUR,8)
TUR.tun             <- indices
TUR.ctrl            <- FLSAM.control(TUR,TUR.tun)

TUR.ctrl@states["catch",]                   <- c(0:5,rep(6,2))
TUR.ctrl@cor.F                              <- 2
TUR.ctrl@catchabilities["SNS",ac(1:6)]      <- c(0,0,1,rep(2,3))        + 101
TUR.ctrl@catchabilities["BTS-ISIS",ac(1:7)] <- c(0,0,1,1,2,2,2)         + 201
TUR.ctrl@catchabilities["NL_LPUE",ac(1)]    <- 0                        + 301
TUR.ctrl@f.vars["catch",]                   <- c(0,0,1,1,2,2,2,2)
TUR.ctrl@logN.vars[]                        <- c(0,rep(1,7))
TUR.ctrl@obs.vars["catch",]                 <- c(0,0,1,1,1,2,2,2)       + 101
TUR.ctrl@obs.vars["SNS",ac(1:6)]            <- c(0,0,1,1,1,1)           + 201
TUR.ctrl@obs.vars["BTS-ISIS",ac(1:7)]       <- c(0,0,0,1,2,2,2)         + 301
TUR.ctrl@obs.vars["NL_LPUE",ac(1)]          <- 0                        + 401
TUR.ctrl@cor.obs[]                          <- NA
TUR.ctrl@cor.obs["SNS",1:5]                 <- c(0,1,1,1,1)
TUR.ctrl@cor.obs.Flag[2]                    <- af("AR")
TUR.ctrl@biomassTreat[4]                    <- 2
TUR.ctrl                                    <- update(TUR.ctrl)
### ------------------------------------------------------------------------------------------------------
###   4. Run assessment
### ------------------------------------------------------------------------------------------------------

TUR.sam             <- FLSAM(TUR,TUR.tun,TUR.ctrl)
TUR.ctrl@residuals  <- FALSE; TUR.sam@control@residuals <- FALSE
TUR.retro           <- retro(TUR,TUR.tun,TUR.ctrl,retro=7,base.assess=TUR.sam)
TUR.ctrl@residuals  <- TRUE; TUR.sam@control@residuals <- TRUE

### ------------------------------------------------------------------------------------------------------
###   5. Diagnostics
### ------------------------------------------------------------------------------------------------------

source(file.path(codePath,"03b_runDiagnostics.r"))
AIC(TUR.sam)
mean(mohns.rho(TUR.retro,ref.year=2016,span=7)$rho[1:7])
mean(mohns.rho(TUR.retro,ref.year=2016,span=7,type="ssb")$rho[1:7])
mean(mohns.rho(TUR.retro,ref.year=2016,span=7,type="rec")$rho[1:7])

#Checks
#default:   AIC = 757.4095  Mohns = -32.5927
#cor.F=2:   AIC = 757.4095  Mohns = -32.5927      Comment: this was default
#cor.F=0:   AIC = 783.149   Mohns = -31.55246
#cor.F=1:   AIC = 781.8074  Mohns = -34.41213     Comment:  6-year peal only as model would not converge with 7
#decision: Cor.F = 2
#q SNS def: AIC = 757.4095  Mohns = -32.5927      Comment: this was default
#q SNS af:  AIC = 761.3539  Mohns = -32.607
#q SNS 01:  AIC = 756.1726  Mohns = -30.72667
#decision: SNS bind = c(0,0,1,2,2,2)
#q BTS def: AIC = 756.1726  Mohns = -30.72667
#q BTS af:  AIC = 758.1933  Mohns = -39.19334
#q BTS 0167 AIC = 754.2291  Mohns = -35.22909
#q BTS 01567AIC = 755.2692  Mohns = -32.2692
#decision: BTS bind = c(0,0,1,1,2,2,2)
#fvars def: AIC =
#fvars af:  AIC = 757.7532  Mohns = -27.11146     Comment: 0 1 2 2 3 3 3 3
#fvars 1:   AIC = 753.8252  Mohns = -29           Comment: 0 0 1 1 2 2 2 2
#fvars 2:   AIC = 749.9858  Mohns = -27.94554     Comment: 0 1 1 1 2 2 2 2
#decision is fvars 1  0 0 1 1 2 2 2 2
#nvars: no testing, nothing worked
#decision to use default c(0,rep(1,7))
#obscatchd: AIC = 754.1348  Mohns = -33.99172     Comment: 0 1 2 2 3 3 3 3
#obscatch1: AIC = 757.1431  Mohns = -29.40655     Comment: 0 0 1 1 1 2 2 2
#obscatch2: AIC = 755.6416  Mohns = -29.50202     Comment: 0 0 1 1 1 2 2 2
#obscatch3: AIC = 756.1726  Mohns = -30.72682     Comment: 0 0 1 1 2 2 3 3
#obscatch4: AIC = 756.5102  Mohns = -32.0598      Comment: 0 1 2 2 2 3 3 3
#decision to use obscatch2
#obssnsd:   AIC = 755.6416  Mohns = -29.50202     Comment: 0 0 1 1 2 2
#obssns1:   AIC = 759.3523  Mohns = -28.67256     Comment: 0 1 2 3 4 4
#obssns2:   AIC = 755.1534  Mohns = -28.54259     Comment: 0 0 1 1 1 1
#decision to use obssns2
#obsbtsd:   AIC = 755.1534  Mohns = -28.54259     Comment: 0 0 1 1 2 2 2
#obsbts1:   AIC = 758.2768  Mohns = -28.02168     Comment: 0 1 2 3 4 5 5
#obsbts2:   AIC = 752.4609  Mohns = -27.43552     Comment: 0 0 0 1 2 2 2
#decision is to use obsbts2
#corsnsd:   AIC = 752.4609  Mohns = -27.43552     Comment: 0 1 1 1 1
#corsnsn:   AIC = 756.6974  Mohns = -29.79415
#corsns1:   AIC = 755.4328  Mohns = -27.79802     Comment: 0 0 1 1 1
#corsns2:   AIC = 756.1605                        Comment: 0 1 2 3 3
#corsns3:   AIC = 754.4207  Mohns = -27.27599     Comment: 0 1 2 2 2
#decision is to use default
