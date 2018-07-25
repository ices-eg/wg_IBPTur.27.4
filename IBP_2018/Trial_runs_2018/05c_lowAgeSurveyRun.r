#-------------------------------------------------------------------------------
# IBPNEW 2012
# Code to generate the .dat file for input into the ADMB model
# Uses Lowestoft format input files
# David Miller, some code adapted from FLR-project
# Date: 3 Oct 2012
#-------------------------------------------------------------------------------
rm(list=ls())

library(FLCore);
library(FLAssess);
library(FLSAM)
library(FLEDA); 

library(mgcv)
library(splines);
library(scales);
library(gplots);
library(grid);
library(gridExtra);
library(latticeExtra)
library(sas7bdat)
library(TMB);

# Set paths to correct folders 
#Path      <- paste0("D:/Bestandsbeheer/WGNSSK/2018/Stock/tur-nsea/10_Interbenchmark/assessment runs/")
Path      <- "D:/Repository/wg_IBPTur.27.4/IBP_2018/"
dataPath  <- paste(Path,"Lowestoft_files/",sep="")
outPath   <- paste(Path,"Trial_runs_2018/output/base_maxAgeSurveys/",sep="")
codePath  <- paste(Path,"Trial_runs_2018/source/",sep="")

assyear <- 2018
retroYr <- 2018

## Source methods/functions
source(paste(codePath,"01_smoothweights.r",sep=""))
source(paste(codePath,"03a_setupStockIndices.r",sep=""))

run       <- "base_maxAgeSurveys"
sens      <- ""

### ------------------------------------------------------------------------------------------------------
###   2. Read and process assessment input data
### ------------------------------------------------------------------------------------------------------

# Not taking IBTS_Q1 in assessment... this was done in the IBP 2017.
names(indices)
# Only select SNS, BTS-ISIS and NL_LPUE_modelD
indices             <- FLIndices(list(window(trim(indices[[1]],age=1:6),start=2004),window(indices[[2]],start=1991),indices[[5]]))
names(indices)      <- c("SNS","BTS-ISIS","NL_LPUE")

### ------------------------------------------------------------------------------------------------------
###   3. Setup data structure for SAM assessment
### ------------------------------------------------------------------------------------------------------

TUR                         <- window(stock,start=1981)
TUR@catch.n[,ac(2000:2002)] <- -1
TUR@landings.n[]            <- TUR@catch.n
TUR.tun                     <- indices
TUR.ctrl                    <- FLSAM.control(TUR,TUR.tun)

### ------------------------------------------------------------------------------------------------------
###   3. Setup data structure for SAM assessment
### ------------------------------------------------------------------------------------------------------

TUR.sams        <- list()
#TUR.sams        <- new("FLSAMs")
TUR.sams.retro  <- list()
for(pg in 7:3){
  TUR@name            <- paste(stock@name,"pgSurvey",pg,sep=" ")
  TUR.tun             <- indices
  for(iTun in c("SNS","BTS-ISIS")){
    if(pg < range(TUR.tun[[iTun]])["max"]){
      TUR.tun[[iTun]]   <- trim(TUR.tun[[iTun]],age=1:pg)
    }
  }
  TUR.ctrl            <- FLSAM.control(TUR,TUR.tun)

  TUR.ctrl@states["catch unique",]              <- c(0:6,rep(7,3))
  TUR.ctrl@cor.F                                <- 2
  if(pg <  range(indices[["SNS"]])["max"])
    TUR.ctrl@catchabilities["SNS",ac(1:pg)]     <- c(0:2,rep(3,3))[1:pg]       + 101
  if(pg >=  range(indices[["SNS"]])["max"])
    TUR.ctrl@catchabilities["SNS",ac(1:6)]      <- c(0:2,rep(3,3))             + 101
  if(pg <  range(indices[["BTS-ISIS"]])["max"])
    TUR.ctrl@catchabilities["BTS-ISIS",ac(1:pg)]<- c(0,0,1,1,rep(2,3))[1:pg]   + 201
  if(pg >=  range(indices[["BTS-ISIS"]])["max"])
    TUR.ctrl@catchabilities["BTS-ISIS",ac(1:7)] <- c(0,0,1,1,rep(2,3))         + 201
  TUR.ctrl@catchabilities["NL_LPUE",ac(1)]      <- 0                           + 301
  #TUR.ctrl@catchabilities["IBTS_Q1",ac(1)]     <- 0                             + 401
  TUR.ctrl@f.vars["catch unique",]              <- c(0,1,2,2,3,3,3,4,4,4)
  TUR.ctrl@logN.vars[]                          <- c(0,rep(1,9))
  TUR.ctrl@obs.vars["catch unique",]            <- c(0,1,2,2,3,3,4,4,4,4)      + 101
  if(pg <  range(indices[["SNS"]])["max"])
   TUR.ctrl@obs.vars["SNS",ac(1:pg)]            <- c(0,0,1,2,3,3)[1:pg]        + 201
  if(pg >=  range(indices[["SNS"]])["max"])
    TUR.ctrl@obs.vars["SNS",ac(1:6)]            <- c(0,0,1,2,3,3)              + 201
  if(pg <  range(indices[["BTS-ISIS"]])["max"])
    TUR.ctrl@obs.vars["BTS-ISIS",ac(1:pg)]      <- c(0,0,0,1,2,3,3)[1:pg]      + 301
  if(pg >=  range(indices[["BTS-ISIS"]])["max"])
    TUR.ctrl@obs.vars["BTS-ISIS",ac(1:7)]       <- c(0,0,0,1,2,3,3)            + 301
  TUR.ctrl@obs.vars["NL_LPUE",ac(1)]            <- 0                           + 401
#  TUR.ctrl@obs.vars["IBTS_Q1",ac(1)]            <- 0                           + 501
  TUR.ctrl@cor.obs[]                            <- NA
  if(pg <  range(indices[["SNS"]])["max"])
    TUR.ctrl@cor.obs["SNS",1:(pg-1)]            <- c(0,rep(1,4))[1:(pg-1)]
  if(pg >=  range(indices[["SNS"]])["max"])
    TUR.ctrl@cor.obs["SNS",1:5]                 <- c(0,rep(1,4))
  TUR.ctrl@cor.obs.Flag[2]                      <- af("AR")
  TUR.ctrl@biomassTreat[4]                      <- 2
  TUR.ctrl                                      <- update(TUR.ctrl)
  ### ------------------------------------------------------------------------------------------------------
  ###   4. Run assessment
  ### ------------------------------------------------------------------------------------------------------

  TUR.sam             <- FLSAM(TUR,TUR.tun,TUR.ctrl)
  TUR.ctrl@residuals  <- FALSE; TUR.sam@control@residuals <- FALSE
  TUR.retro           <- retro(TUR,TUR.tun,TUR.ctrl,retro=5) #base.assess=TUR.sam
  TUR.sams[[ac(pg)]]      <- TUR.sam
  TUR.sams.retro[[ac(pg)]]<- TUR.retro
}

TUR.sams <- as(TUR.sams,"FLSAMs")
plot(TUR.sams)
TUR.sams.retro <- as(TUR.sams.retro,"FLSAMs")
plot(TUR.sams.retro)

### ------------------------------------------------------------------------------------------------------
###   5. Diagnostics
### ------------------------------------------------------------------------------------------------------
#- Mohns rho
lapply(lapply(TUR.sams.retro,mohns.rho,ref.year=2017,span=5,type="ssb"),function(x){return(mean(x$rho[1:5]))})

for(pg in 7:3){
  sens <- paste0("pgSurvey",pg)
  TUR.sam <- TUR.sams[[ac(pg)]]
  TUR.retro <- TUR.sams.retro[[ac(pg)]]
  source(file.path(codePath,"03b_runDiagnostics.r"))
}

pdf(file.path(outPath,paste0(run,"_","maxAgeSurveyscomb","assessmentOut.pdf")))
plot(TUR.sams)

par(mfrow=c(2,1))
print(plot(AIC(TUR.sams),ylab="AIC",xaxt="n",las=2,pch=19,xlab=""))
axis(1,at=1:5,labels=paste0("pg ",names(TUR.sams)),las=1)
print(grid())

print(plot(unlist(lapply(lapply(TUR.sams.retro,mohns.rho,ref.year=2016,span=5,type="fbar"),function(x){return(mean(x$rho[1:5]))})),xlab="",ylab="Mohns rho (5-year peel)",xaxt="n",las=2,pch=19))
axis(1,at=1:5,labels=paste0("pg ",names(TUR.sams.retro)))
print(grid())
dev.off()

