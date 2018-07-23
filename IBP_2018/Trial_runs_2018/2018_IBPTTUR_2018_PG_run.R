#-------------------------------------------------------------------------------
# WGNSSK 2018 prep - Final run JB
# Uses Lowestoft 2017 format input files
# Changes made to code using spline to estimate weca and west. 
#-------------------------------------------------------------------------------
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
#Path      <- "D:/Repository/Turbot/IBP_2018/"

# dataPath  <- paste(Path,"Lowestoft_files/",sep="")
# outPath   <- paste(Path,"Trial_runs_2018/output/",sep="")
# codePath  <- paste(Path,"Trial_runs_2018/source/",sep="")

# to check whether there is something in the source codes wrong - shouldnt be as they are copies
# Path      <- paste0("D:/Bestandsbeheer/WGNSSK/2018/Stock/tur-nsea/")
# dataPath  <- paste(Path,"4_Lowestoft_files/",sep="")
# outPath   <- "D:/Bestandsbeheer/WGNSSK/2018/Stock/tur-nsea/10_Interbenchmark/assessment runs/5_assessment/output/"
# codePath  <- paste(Path,"5_assessment/source/",sep="")

# assyear <- 2017
# stock_Name  <- "tur-nsea"
run       <- "pgroup"
sens      <- ""

## Assessment settings
# Year (= year when assessment conducted.  i.e. have data up to assYear-1)
# assyear <- 2018
# retroYr <- 2018
# last year of data:
# endYear <- assyear-1  

## Source methods/functions
source(paste(codePath,"01_smoothweights.r",sep=""))
source(paste(codePath,"03a_setupStockIndices.r",sep=""))

### ------------------------------------------------------------------------------------------------------
###   2. Read and process assessment input data
### ------------------------------------------------------------------------------------------------------

names(indices)
# Only select SNS, BTS-ISIS and NL_LPUE_modelD
indices             <- FLIndices(list(window(trim(indices[[1]],age=1:6),start=2004),window(indices[[2]],start=1991),indices[[5]]))
names(indices)      <- c("SNS","BTS-ISIS","NL_LPUE")

### ------------------------------------------------------------------------------------------------------
###   3. Setup data structure for SAM assessment
### ------------------------------------------------------------------------------------------------------

TUR.sams  <- new("FLSAMs")
TUR.sams.retro <- list()
for(pg in 7:6){
  TUR                         <- window(stock,start=1981)
  TUR@catch.n[,ac(2000:2002)] <- -1
  TUR@landings.n[]            <- TUR@catch.n
  TUR                         <- setPlusGroup(TUR,pg)
  TUR@name                    <- paste(TUR@name,"pg",pg,sep=" ")
  TUR.tun                     <- indices
  if(pg == 6){
    index                 <- TUR.tun[["BTS-ISIS"]]@index
    index[pg,]            <- quantSums(index[pg:7,])
    TUR.tun[["BTS-ISIS"]] <- trim(TUR.tun[["BTS-ISIS"]],age=1:6)
    TUR.tun[["BTS-ISIS"]]@index[pg,] <- index[pg,]
  }
  TUR.ctrl                    <- FLSAM.control(TUR,TUR.tun)
  
  TUR.ctrl@states["catch unique",]            <- c(0:(pg-2),(pg-2))
  TUR.ctrl@cor.F                              <- 2
  TUR.ctrl@catchabilities["SNS",ac(1:6)]      <- c(0:2,rep(3,3))          + 101
  if(pg > 6)
    TUR.ctrl@catchabilities["BTS-ISIS",ac(1:7)] <- c(0,0,1,1,2,2,2)       + 201
  if(pg == 6)
    TUR.ctrl@catchabilities["BTS-ISIS",ac(1:6)] <- c(0,0,1,1,2,2)         + 201
  TUR.ctrl@catchabilities["NL_LPUE",ac(1)]    <- 0                        + 301
  TUR.ctrl@f.vars["catch unique",]            <- c(0,1,2,2,3,3,4,4,4,4)[1:pg]
  TUR.ctrl@f.vars["catch unique",ncol(TUR.ctrl@f.vars)] <-   TUR.ctrl@f.vars["catch unique",ncol(TUR.ctrl@f.vars)-1]
  TUR.ctrl@logN.vars[]                        <- c(0,rep(1,9))[1:pg]
  TUR.ctrl@obs.vars["catch unique",]          <- c(0,1,2,2,3,3,4,4,4,4)[1:pg]     + 101
  TUR.ctrl@obs.vars["SNS",ac(1:6)]            <- c(0,0,1,2,3,3)         + 201
  if(pg > 6)
    TUR.ctrl@obs.vars["BTS-ISIS",ac(1:7)]     <- c(0,0,0,1,2,3,3)       + 301
  if(pg == 6)
    TUR.ctrl@obs.vars["BTS-ISIS",ac(1:6)]     <- c(0,0,0,1,2,3)         + 301
  TUR.ctrl@obs.vars["NL_LPUE",ac(1)]          <- 0                      + 401
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
  TUR.retro           <- retro(TUR,TUR.tun,TUR.ctrl,retro=5) # base.assess=TUR.sam
  TUR.sams[[ac(pg)]]      <- TUR.sam
  TUR.sams.retro[[ac(pg)]]<- TUR.retro
}

### ------------------------------------------------------------------------------------------------------
###   5. Diagnostics
### ------------------------------------------------------------------------------------------------------
#- Mohns rho
lapply(lapply(TUR.sams.retro,mohns.rho,ref.year=2017,span=6,type="fbar"),function(x){return(mean(x$rho[1:6]))})
lapply(lapply(TUR.sams.retro,mohns.rho,ref.year=2017,span=6,type="ssb"),function(x){return(mean(x$rho[1:6]))})
lapply(lapply(TUR.sams.retro,mohns.rho,ref.year=2017,span=6,type="rec"),function(x){return(mean(x$rho[1:6]))})

for(pg in 7:6){
  sens <- paste0("pg",pg)
  TUR.sam <- TUR.sams[[ac(pg)]]
  TUR.retro <- TUR.sams.retro[[ac(pg)]]
  source(paste0(codePath,"03b_runDiagnostics.r"))
}

pdf(paste0(outPath,"Output",paste0(run,"_","pgcomb","assessmentOut.pdf")))
plot(TUR.sams)
dev.off()
