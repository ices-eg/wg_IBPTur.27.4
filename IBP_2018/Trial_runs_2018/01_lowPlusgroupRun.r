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
Path      <- "W:\\IMARES\\data\\ICES-WG\\IBPTURBOT\\2017\\assessment runs\\"
dataPath  <- paste(Path,"Lowestoft files\\",sep="")
outPath   <- paste(Path,"trial_runs_2017\\",sep="")

## Source methods/functions
source(paste(Path,"nsea_functions.r",sep=""))
source(paste(outPath,"03a_setupStockIndices.r",sep=""))

run       <- "pgroup"
sens      <- ""

### ------------------------------------------------------------------------------------------------------
###   2. Read and process assessment input data
### ------------------------------------------------------------------------------------------------------

indices             <- FLIndices(list(window(trim(indices[[1]],age=1:6),start=2004),window(indices[[2]],start=1991),indices[[3]]))

### ------------------------------------------------------------------------------------------------------
###   3. Setup data structure for SAM assessment
### ------------------------------------------------------------------------------------------------------

TUR.sams  <- new("FLSAMs")
TUR.sams.retro <- list()
for(pg in 7:6){
  TUR                 <- setPlusGroup(stock,pg)
  TUR@name            <- paste(TUR@name,"pg",pg,sep=" ")
  TUR.tun             <- indices
  if(pg == 6){
    index             <- TUR.tun[["BTS-ISIS"]]@index
    index[pg,]        <- quantSums(index[pg:7,])
    TUR.tun[["BTS-ISIS"]] <- trim(TUR.tun[["BTS-ISIS"]],age=1:6)
    TUR.tun[["BTS-ISIS"]]@index[pg,] <- index[pg,]
  }
  TUR.ctrl            <- FLSAM.control(TUR,TUR.tun)

  TUR.ctrl@states["catch unique",]                   <- c(0:(pg-2),(pg-2))
  TUR.ctrl@cor.F                              <- 2
  TUR.ctrl@catchabilities["SNS",ac(1:6)]      <- c(0:2,rep(3,3))          + 101
  if(pg > 6)
    TUR.ctrl@catchabilities["BTS-ISIS",ac(1:7)]<- c(0,0,1,1,rep(2,3))     + 201
  if(pg == 6)
    TUR.ctrl@catchabilities["BTS-ISIS",ac(1:6)]<- c(0,0,1,1,rep(2,2))     + 201
  TUR.ctrl@catchabilities["NL_LPUE",ac(1)]    <- 0                        + 301
  TUR.ctrl@f.vars["catch unique",]                   <- c(0,1,2,2,3,3,3,4,4,4)[1:pg]
  TUR.ctrl@f.vars["catch unique",ncol(TUR.ctrl@f.vars)] <-   TUR.ctrl@f.vars["catch unique",ncol(TUR.ctrl@f.vars)-1]
  TUR.ctrl@logN.vars[]                        <- c(0,rep(1,9))[1:pg]
  TUR.ctrl@obs.vars["catch unique",]                 <- c(0,1,2,2,3,3,4,4,4,4)[1:pg]   + 101
  TUR.ctrl@obs.vars["SNS",ac(1:6)]            <- c(0,0,1,2,3,3)           + 201
  if(pg > 6)
    TUR.ctrl@obs.vars["BTS-ISIS",ac(1:7)]     <- c(0,0,0,1,2,3,3)         + 301
  if(pg == 6)
    TUR.ctrl@obs.vars["BTS-ISIS",ac(1:6)]     <- c(0,0,0,1,2,3)           + 301
  TUR.ctrl@obs.vars["NL_LPUE",ac(1)]          <- 0                        + 401
  TUR.ctrl@cor.obs[]                          <- NA
  TUR.ctrl@cor.obs["SNS",1:5]                 <- c(0,rep(1,4))
  TUR.ctrl@cor.obs.Flag[2]                    <- af("AR")
  TUR.ctrl@biomassTreat[4]                    <- 2
  TUR.ctrl                                    <- update(TUR.ctrl)
  ### ------------------------------------------------------------------------------------------------------
  ###   4. Run assessment
  ### ------------------------------------------------------------------------------------------------------

  TUR.sam             <- FLSAM(TUR,TUR.tun,TUR.ctrl)
  TUR.ctrl@residuals  <- FALSE; TUR.sam@control@residuals <- FALSE
  TUR.retro           <- retro(TUR,TUR.tun,TUR.ctrl,retro=7,base.assess=TUR.sam)
  TUR.sams[[ac(pg)]]      <- TUR.sam
  TUR.sams.retro[[ac(pg)]]<- TUR.retro
}
### ------------------------------------------------------------------------------------------------------
###   5. Diagnostics
### ------------------------------------------------------------------------------------------------------
#- Mohns rho
lapply(lapply(TUR.sams.retro,mohns.rho,ref.year=2016,span=7),function(x){return(mean(x$rho[1:7]))})
for(pg in 10:6){
  sens <- paste0("pg",pg)
  TUR.sam <- TUR.sams[[ac(pg)]]
  TUR.retro <- TUR.sams.retro[[ac(pg)]]
  source(file.path(outPath,"03b_runDiagnostics.r"))
}

pdf(file.path(outPath,"Output",paste0(run,"_","pgcomb","assessmentOut.pdf")))
plot(TUR.sams)
dev.off()

