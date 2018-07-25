#-------------------------------------------------------------------------------
# IBPNEW 2012
# Code to generate the .dat file for input into the ADMB model
# Uses Lowestoft format input files
# David Miller, some code adapted from FLR-project
# Date: 3 Oct 2012
rm(list=ls())

# FLR
# install.packages(pkgs="FLAssess",repos="http://flr-project.org/R")
# install.packages(pkgs="FLEDA",repos="http://flr-project.org/R")
#  install.packages(pkgs="FLCore",repos="http://flr-project.org/R")
# devtools::install_github("fishfollower/SAM/stockassessment", ref="components")
# devtools::install_github("flr/FLSAM", ref="develop_V2")

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

# Set paths to folders (10_Interbenchmark)
#Path      <- paste0("D:/Bestandsbeheer/WGNSSK/2018/Stock/tur-nsea/10_Interbenchmark/assessment runs/")
Path      <- "D:/Repository/wg_IBPTur.27.4/IBP_2018/"
Path      <- "D:/Repository/Turbot/IBP_2018/"

dataPath  <- paste(Path,"Lowestoft_files/",sep="")
outPath   <- paste(Path,"Trial_runs_2018/output/baserun/",sep="")
codePath  <- paste(Path,"Trial_runs_2018/source/",sep="")

## Source methods/functions
source(paste(codePath,"01_smoothweights.r",sep=""))
source(paste(codePath,"03a_setupStockIndices.r",sep=""))

run       <- "data_base"
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

TUR                 <- window(stock,start=1981)
TUR@catch.n[,ac(2000:2002)] <- -1; TUR@landings.n[] <- TUR@catch.n
TUR.tun             <- indices
TUR.ctrl            <- FLSAM.control(TUR,TUR.tun)

TUR.ctrl@states["catch unique",]                   <- c(0:6,rep(7,3))
TUR.ctrl@cor.F                              <- 2
TUR.ctrl@catchabilities["SNS",ac(1:6)]      <- c(0:2,rep(3,3))          + 101
TUR.ctrl@catchabilities["BTS-ISIS",ac(1:7)] <- c(0,0,1,1,rep(2,3))      + 201
TUR.ctrl@catchabilities["NL_LPUE",ac(1)]    <- 0                        + 301
# TUR.ctrl@catchabilities["IBTS_Q1",ac(1)]    <- 0                        + 401
TUR.ctrl@f.vars["catch unique",]                   <- c(0,1,2,2,3,3,3,4,4,4)
TUR.ctrl@logN.vars[]                        <- c(0,rep(1,9))
TUR.ctrl@obs.vars["catch unique",]                 <- c(0,0,1,1,2,2,3,3,4,4)   + 101
TUR.ctrl@obs.vars["SNS",ac(1:6)]            <- c(0,0,1,1,2,2)           + 201
TUR.ctrl@obs.vars["BTS-ISIS",ac(1:7)]       <- c(0,0,1,1,2,2,2)         + 301
TUR.ctrl@obs.vars["NL_LPUE",ac(1)]          <- 0                        + 401
# TUR.ctrl@obs.vars["IBTS_Q1",ac(1)]          <- 0                        + 501
#TUR.ctrl@cor.obs[]                          <- NA
#TUR.ctrl@cor.obs["SNS",1:5]                 <- c(0,rep(1,4))
#TUR.ctrl@cor.obs.Flag[2]                    <- af("AR")
TUR.ctrl@biomassTreat[4]                    <- 2
TUR.ctrl                                    <- update(TUR.ctrl)

### ------------------------------------------------------------------------------------------------------
###   4. Run assessment
### ------------------------------------------------------------------------------------------------------

TUR.sam             <- FLSAM(TUR,TUR.tun,TUR.ctrl)
TUR.ctrl@residuals  <- FALSE; TUR.sam@control@residuals <- FALSE
TUR.retro           <- retro(TUR,TUR.tun,TUR.ctrl,retro=5)
source(file.path(codePath,"03b_runDiagnostics.r"))

#- LOA
TUR.sams   <- list()
TUR.retros <- list()
for(i in names(TUR.tun)){
  TUR.ctrlLOA         <- drop.from.control(TUR.ctrl,fleets=i)
  TUR.LOA             <- TUR
  idxkeep             <- which(names(TUR.tun)!=i)
  TUR.tunLOA          <- TUR.tun[idxkeep]
  TUR.ctrlLOA@residuals  <- TRUE
  TUR.sams[[i]]       <- FLSAM(TUR.LOA,TUR.tunLOA,TUR.ctrlLOA)
  TUR.ctrlLOA@residuals  <- FALSE
  TUR.retros[[i]]     <- retro(TUR.LOA,TUR.tunLOA,TUR.ctrlLOA,retro=5)
}

TUR.sams   <- as(TUR.sams,"FLSAMs")
plot(TUR.sams)

TUR.retros <- as(TUR.retros,"FLSAMs")
plot(TUR.retros[[1]])

### ------------------------------------------------------------------------------------------------------
###   5. Diagnostics
### ------------------------------------------------------------------------------------------------------


for(i in names(TUR.tun)){
  sens <- paste0(run,"_LOA_",i)
  TUR.sam <- TUR.sams[[i]]
  TUR.retros <- TUR.retros[[i]]
  source(file.path(codePath,"03b_runDiagnostics.r"))
}


pdf(paste0(outPath,paste0(run,"_","LOA_","comb","assessmentOut.pdf")))
plot(TUR.sams)

par(mfrow=c(2,1))
print(plot(AIC(TUR.sams),ylab="AIC",xaxt="n",las=2,pch=19,xlab=""))
axis(1,at=1:3,labels=names(TUR.sams),las=1)
print(grid())

print(plot(unlist(lapply(lapply(TUR.retros,mohns.rho,ref.year=2017,span=5,type="ssb"),function(x){return(mean(x$rho[1:5]))})),xlab="",ylab="Mohns rho (5-year peel)",xaxt="n",las=2,pch=19,main="SSB"))
axis(1,at=1:3,labels=paste0("pg ",names(TUR.sams.retro)))
print(grid())

print(plot(unlist(lapply(lapply(TUR.sams.retro,mohns.rho,ref.year=2017,span=5,type="fbar"),function(x){return(mean(x$rho[1:5]))})),xlab="",ylab="Mohns rho (5-year peel)",xaxt="n",las=2,pch=19,main="Fbar"))
axis(1,at=1:3,labels=paste0("pg ",names(TUR.sams.retro)))
print(grid())

print(plot(unlist(lapply(lapply(TUR.sams.retro,mohns.rho,ref.year=2017,span=5,type="rec"),function(x){return(mean(x$rho[1:5]))})),xlab="",ylab="Mohns rho (5-year peel)",xaxt="n",las=2,pch=19,main="Recruitment"))
axis(1,at=1:3,labels=paste0("pg ",names(TUR.sams.retro)))
print(grid())
dev.off()

