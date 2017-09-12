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

run       <- "M"
sens      <- ""

### ------------------------------------------------------------------------------------------------------
###   2. Read and process assessment input data
### ------------------------------------------------------------------------------------------------------

indices             <- FLIndices(list(window(trim(indices[[1]],age=1:6),start=2004),window(indices[[2]],start=1991),indices[[3]]))

### ------------------------------------------------------------------------------------------------------
###   3. Setup data structure for SAM assessment
### ------------------------------------------------------------------------------------------------------

TUR                 <- stock
TUR.tun             <- indices
TUR.ctrl            <- FLSAM.control(TUR,TUR.tun)

TUR.ctrl@states["catch",]                   <- c(0:6,rep(7,3))
TUR.ctrl@cor.F                              <- 2
TUR.ctrl@catchabilities["SNS",ac(1:6)]      <- c(0:2,rep(3,3))          + 101
TUR.ctrl@catchabilities["BTS-ISIS",ac(1:7)] <- c(0,0,1,1,rep(2,3))      + 201
TUR.ctrl@catchabilities["NL_LPUE",ac(1)]    <- 0                        + 301
TUR.ctrl@f.vars["catch",]                   <- c(0,1,2,2,3,3,3,4,4,4)
TUR.ctrl@logN.vars[]                        <- c(0,rep(1,9))
TUR.ctrl@obs.vars["catch",]                 <- c(0,1,2,2,3,3,4,4,4,4)   + 101
TUR.ctrl@obs.vars["SNS",ac(1:6)]            <- c(0,0,1,2,3,3)           + 201
TUR.ctrl@obs.vars["BTS-ISIS",ac(1:7)]       <- c(0,0,0,1,2,3,3)         + 301
TUR.ctrl@obs.vars["NL_LPUE",ac(1)]          <- 0                        + 401
TUR.ctrl@cor.obs[]                          <- NA
TUR.ctrl@cor.obs["SNS",1:5]                 <- c(0,rep(1,4))
TUR.ctrl@cor.obs.Flag[2]                    <- af("AR")
TUR.ctrl@biomassTreat[4]                    <- 2
TUR.ctrl                                    <- update(TUR.ctrl)
### ------------------------------------------------------------------------------------------------------
###   4. Run assessment
### ------------------------------------------------------------------------------------------------------

TUR.sams  <- new("FLSAMs")
TUR.sams.retro <- list()
for(iM in c(0.2*0.75,0.2*(1/0.75),"Lorenzen","LorenzenS")){
  name    <- paste("M",iM)
  if(!iM %in% c("Lorenzen","LorenzenS"))
    TUR@m[]         <- an(iM)
  if(iM == "Lorenzen")
    TUR@m[]         <- 3.69*(TUR@stock.wt*1000)^-0.305
  if(iM == "LorenzenS")
    TUR@m[]         <- sweep(3.69*(TUR@stock.wt*1000)^-0.305,2:6,(quantMeans(3.69*(TUR@stock.wt*1000)^-0.305)/0.2),"/")

  TUR.sam             <- FLSAM(TUR,TUR.tun,TUR.ctrl)
  TUR.ctrl@residuals  <- FALSE; TUR.sam@control@residuals <- FALSE
  TUR.retro           <- retro(TUR,TUR.tun,TUR.ctrl,retro=7,base.assess=TUR.sam)
  TUR.sams[[name]]    <- TUR.sam
  TUR.sams.retro[[name]] <- TUR.retro
  TUR.ctrl@residuals  <- TRUE; TUR.sam@control@residuals <- TRUE
}


### ------------------------------------------------------------------------------------------------------
###   5. Diagnostics
### ------------------------------------------------------------------------------------------------------
for(iM in c(0.2*0.75,0.2*(1/0.75),"Lorenzen","LorenzenS")){
  sens <- paste0("M_",iM)
  name    <- paste("M",iM)
  TUR.sam <- TUR.sams[[name]]
  TUR.retro <- TUR.sams.retro[[name]]
  source(file.path(codePath,"03b_runDiagnostics.r"))
}


pdf(file.path(outPath,paste0(run,"_","Mcomb","assessmentOut.pdf")))
plot(TUR.sams)

par(mfrow=c(2,1))
print(plot(AIC(TUR.sams),ylab="AIC",xaxt="n",las=2,pch=19,xlab=""))
axis(1,at=1:4,labels=names(TUR.sams),las=1)
print(grid())
print(plot(unlist(lapply(lapply(TUR.sams.retro,mohns.rho,ref.year=2016,span=7),function(x){return(mean(x$rho[1:7]))})),xlab="",ylab="Mohns rho (7-year peel)",xaxt="n",las=2,pch=19))
axis(1,at=1:4,labels=names(TUR.sams.retro))
print(grid())
dev.off()

