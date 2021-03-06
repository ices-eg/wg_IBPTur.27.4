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
Path      <- "D:/Repository/wg_IBPTur.27.4/IBP_2018/"
#Path      <- "D:/Repository/Turbot/IBP_2018/"

dataPath  <- paste(Path,"Lowestoft_files/",sep="")
outPath   <- paste(Path,"Trial_runs_2018/output/parsettings",sep="")
codePath  <- paste(Path,"Trial_runs_2018/source/",sep="")

## Source methods/functions
source(paste(codePath,"01_smoothweights.r",sep=""))
source(paste(codePath,"03a_setupStockIndices.r",sep=""))

run       <- "final_"
sens      <- "IBP_2018"

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

TUR                         <- window(stock,start=1981)
TUR@catch.n[,ac(2000:2002)] <- -1
TUR@landings.n[]            <- TUR@catch.n
TUR                         <- setPlusGroup(TUR,8)
TUR.tun                     <- indices
TUR.ctrl                    <- FLSAM.control(TUR,TUR.tun)

TUR.ctrl@states["catch unique",]            <- c(0,1,2,3,4,5,6,6)
TUR.ctrl@cor.F                              <- 2
TUR.ctrl@catchabilities["SNS",ac(1:6)]      <- c(0,0,1,2,2,2)           + 101
TUR.ctrl@catchabilities["BTS-ISIS",ac(1:7)] <- c(0,0,1,1,2,2,2)         + 201
TUR.ctrl@catchabilities["NL_LPUE",ac(1)]    <- 0                        + 301
TUR.ctrl@f.vars["catch unique",]            <- c(0,1,2,2,3,3,4,4)
TUR.ctrl@logN.vars[]                        <- c(0,rep(1,7))
TUR.ctrl@obs.vars["catch unique",]          <- c(0,1,2,2,3,3,4,4)       + 101
TUR.ctrl@obs.vars["SNS",ac(1:6)]            <- c(0,0,1,2,2,2)           + 201
TUR.ctrl@obs.vars["BTS-ISIS",ac(1:7)]       <- c(0,0,0,1,2,2,2)         + 301
TUR.ctrl@obs.vars["NL_LPUE",ac(1)]          <- 0                        + 401
TUR.ctrl@cor.obs[]                          <- NA
#TUR.ctrl@cor.obs["BTS-ISIS",1:6]           <- c(0,1,1,1,1,1)
TUR.ctrl@cor.obs["SNS",1:5]                <- c(0,0,0,0,0)
TUR.ctrl@cor.obs.Flag[2]                   <- af("AR")
TUR.ctrl@biomassTreat[4]                    <- 2
#TUR.ctrl@residuals                          <- FALSE
TUR.ctrl                                    <- update(TUR.ctrl)

### ------------------------------------------------------------------------------------------------------
###   4. Run assessment
### ------------------------------------------------------------------------------------------------------

TUR.sam             <- FLSAM(TUR,TUR.tun,TUR.ctrl)
TUR.ctrl@residuals  <- FALSE; TUR.sam@control@residuals <- FALSE
TUR.retro           <- retro(TUR,TUR.tun,TUR.ctrl,retro=5)
source(file.path(codePath,"03b_runDiagnostics.r"))

plot(TUR.sam)
plot(TUR.retro)
#save(TUR.sam, file = paste0(outPath,"TUR.sam_2018_assesment.Rdata", sep = ""))

#- Manual retro for retro in stockweight
load(file.path(outPath,"sswtsRetro.RData"))
load(file.path(outPath,"scwtsRetro.RData"))
TUR.retroWeight     <- list()
TUR.retroWeight[[ac(2017)]] <- TUR.sam
for(iYr in 2016:2012){
  TUR.tun.temp <- window(TUR.tun,end=iYr)
  TUR.ctrl.temp <- TUR.ctrl
  TUR                         <- window(stock,start=1981,end=iYr)
  TUR@catch.n[,ac(2000:2002)] <- -1
  TUR@landings.n[]            <- TUR@catch.n
  TUR@catch.wt[]              <- get(paste0("scwts",iYr))[,-c(1:6)] #trim off years 1975:1981
  TUR@stock.wt[]              <- get(paste0("sswts",iYr))[,-c(1:6)] #trim off years 1975:1981
  TUR                         <- setPlusGroup(TUR,9)
  TUR.temp <- TUR
  TUR.ctrl.temp@name <- as.character(iYr)
  TUR.ctrl.temp@range["maxyear"] <- max(TUR.temp@range["maxyear"],
                                    max(sapply(TUR.tun.temp,function(x) max(x@range[c("maxyear")]))))
  TUR.retroWeight[[ac(iYr)]] <- FLSAM(TUR.temp,TUR.tun.temp,TUR.ctrl.temp)
}
TUR.retroWeight <- as(TUR.retroWeight,"FLSAMs")
mean(mohns.rho(TUR.retroWeight,ref.year=2017,span=5,type="ssb")$rho[1:5])
mean(mohns.rho(TUR.retroWeight,ref.year=2017,span=5,type="fbar")$rho[1:5])
mean(mohns.rho(TUR.retroWeight,ref.year=2017,span=5,type="rec")$rho[1:5])
save(TUR.retroWeight,file=file.path(outPath,"TUR.retroWeight.RData"))


### ------------------------------------------------------------------------------------------------------
###   5. Diagnostics
### ------------------------------------------------------------------------------------------------------

obsvar.plot(TUR.sam)


# catch <- catchabilities(TUR.sam)
# print(xyplot(value+ubnd+lbnd ~ age | fleet,catch,
#              scale=list(alternating=FALSE,y=list(relation="free")),as.table=TRUE,
#              type="l",lwd=c(2,1,1),col=c("black","grey","grey"),
#              subset=fleet %in% c("SNS","BTS-ISIS","NL_LPUE_age"),
#              main="Survey catchability parameters",ylab="Catchability",xlab="Age"))

AIC(TUR.sam)

# mohns rho
mean(mohns.rho(TUR.retro,ref.year=2017,span=5,type="ssb")$rho[1:5])
mean(mohns.rho(TUR.retro,ref.year=2017,span=5,type="fbar")$rho[1:5])
mean(mohns.rho(TUR.retro,ref.year=2017,span=5,type="rec")$rho[1:5])

#Checks (zie excel)



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