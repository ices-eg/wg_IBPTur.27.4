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

# Set paths to folders
Path      <- "W:\\IMARES\\data\\ICES-WG\\IBPTURBOT\\2017\\assessment runs\\"
dataPath  <- paste(Path,"Lowestoft files\\",sep="")
outPath   <- paste(Path,"trial_runs_2017\\",sep="")

## Source methods/functions
source(paste(Path,"nsea_functions.r",sep=""))

### ------------------------------------------------------------------------------------------------------
###  1. Settings
### ------------------------------------------------------------------------------------------------------

## Assessment settings
## Stock name
stock_Name          <- "tur-nsea"
# Year (= year when assessment conducted.  i.e. have data up to assYear-1)
assYear             <- 2017
retroYr             <- 2017
endYear             <- min(assYear,retroYr)-1
Units               <- c("tonnes","thousands","kg")     # totals (product of the next two), numbers, weights-at-age
maxA                <- 10
pGrp                <- T
minFbar             <- 2
maxFbar             <- 6


#------ RUN
run                 <- "LPUE_UK_add"
#------ SENS
sens                <- "trim89_corF2"

### ------------------------------------------------------------------------------------------------------
###   2. Read and process assessment input data
### ------------------------------------------------------------------------------------------------------

## Read stock data
stock               <- readFLStock(paste(dataPath,"index_raw.txt", sep=""))
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
indices             <- readFLIndices(paste(dataPath,"fleet_tmb_longer_LPUE.txt", sep=""), na.strings="-1")
if(run == "base")
  indices           <- FLIndices(list(window(trim(indices[[1]],age=1:6),start=2004),window(indices[[2]],start=1991),indices[[3]]))
if(run == "trim89")
  indices           <- FLIndices(list(window(trim(indices[[1]],age=1:6),start=2004),window(indices[[2]],start=1991),indices[[3]]))
if(run == "trim98")
  indices           <- FLIndices(list(window(trim(indices[[1]],age=1:6),start=2004),window(indices[[2]],start=1998),window(indices[[3]],start=1998)))
if(run == "LPUEmin6"){
  trimLPUE          <- quantSums(trim(indices[[4]],age=3:9)@index)
  indices           <- FLIndices(list(window(trim(indices[[1]],age=1:6),start=2004),window(indices[[2]],start=1998),indices[[3]]))
  indices[[3]]@index<- trimLPUE}
if(run == "LPUE_UK_add")
  indices           <- FLIndices(list(window(trim(indices[[1]],age=1:6),start=2004),window(indices[[2]],start=1991),indices[[3]],window(indices[[5]],start=1991)))
if(run == "LPUE_UK_alone")
  indices           <- FLIndices(list(window(trim(indices[[1]],age=1:6),start=2004),window(indices[[2]],start=1991),window(indices[[5]],start=1991)))



indexVals           <- lapply(indices, index)
numIndices          <- length(indexVals)
indMPs              <- list()
for(ind in names(indices))
  indMPs[[ind]]     <- as.numeric((range(indices[[ind]])["startf"]+range(indices[[ind]])["endf"])/2)

#- Update catch and stock weights based on spline fitted model
load(file.path(outPath,"smoothedWeights.RData"))
stock@catch.wt[]    <- scwts
stock@landings.wt[] <- stock@catch.wt
stock@discards.wt[] <- stock@catch.wt
stock@stock.wt[]    <- sswts

# Correct for SOP in landings.n
landings.n(stock)@.Data[which(landings.n(stock)==-1)] <- NA
landings.n(stock)   <- sweep(landings.n(stock),2,computeLandings(stock)/landings(stock),"/")
landings.n(stock)@.Data[which(is.na(landings.n(stock)))] <- -1
catch.n(stock)      <- landings.n(stock)

if(run == "trim89" | sens == "trim89")
  stock             <- window(stock,start=1989)
if(run == "trim98" | sens == "trim98")
  stock             <- window(stock,start=1998)

### ------------------------------------------------------------------------------------------------------
###   3. Setup data structure for SAM assessment
### ------------------------------------------------------------------------------------------------------

surveys             <- lapply(indices,index)
surveys             <- lapply(surveys,function(x){x[,drop=T]})
surveys             <- lapply(surveys,function(x){if(class(x)=="matrix"){ret <- t(x)} else { ret <- t(t(x))}; return(ret)})
for(i in 1:length(surveys))
  attr(surveys[[i]],"time") <- lapply(indices,function(x){return(range(x)[c("startf","endf")])})[[i]]
if(run=="LPUE_UK_add")
  names(surveys)      <- c("SNS","BTS","LPUE_NL","LPUE_UK")
if(run=="LPUE_UK_alone")
  names(surveys)      <- c("SNS","BTS","LPUE_UK")
if(length(grep("UK",run))==0)
  names(surveys)      <- c("SNS","BTS","LPUE_NL")
for(i in grep("LPUE",names(surveys)))
  colnames(surveys[[i]]) <- -1

#- Combine survey, catch and population information into dat file
dat                 <- setup.sam.data(surveys=surveys,
                                      residual.fleet=   t(stock@catch.n[,drop=T]),
                                      prop.mature=      t(stock@mat[,drop=T]),
                                      stock.mean.weight=t(stock@stock.wt[,drop=T]),
                                      catch.mean.weight=t(stock@catch.wt[,drop=T]),
                                      dis.mean.weight=  t(stock@catch.wt[,drop=T]),
                                      land.mean.weight= t(stock@landings.wt[,drop=T]),
                                      prop.f=           t(stock@harvest.spwn[,drop=T]),
                                      prop.m=           t(stock@m.spwn[,drop=T]),
                                      natural.mortality=t(stock@m[,drop=T]),
                                      land.frac=        t((stock@landings.n/stock@catch.n)[,drop=T]))

setminonetoNA       <- function(x){for(i in 1:length(x)){idx <- which(x[[i]]==-1);x[[i]][idx] <- NA};return(x)}
setNAtominone       <- function(x){for(i in 1:length(x)){idx <- which(is.na(x[[i]]));x[[i]][idx] <- -1};return(x)}

#- Generate default configuration file
conf                <- defcon(dat)
conf                <- setminonetoNA(conf)
#- Set the number of F random walks (estimated Fs at age)
conf$keyLogFsta[1,] <- c(0:6,rep(7,3))
#- Set the random walk correlation flag
conf$corFlag[]      <- 2
#- Set the catchability parameter bindings
conf$keyLogFpar[2,] <- c(0:2,rep(3,3),rep(NA,4))
conf$keyLogFpar[3,] <- c(0,0,1,1,rep(2,3),rep(NA,3))  + max(conf$keyLogFpar[2,],na.rm=T) + 1
conf$keyLogFpar[4,] <- c(0,rep(NA,9))                 + max(conf$keyLogFpar[3,],na.rm=T) + 1
if(run == "LPUE_UK_add")
  conf$keyLogFpar[5,] <- c(0,rep(NA,9))               + max(conf$keyLogFpar[4,],na.rm=T) + 1
#- Set random walk F variances
conf$keyVarF[1,]    <- c(0,1,2,2,3,3,3,4,4,4)
#- Set random walk N variances
conf$keyVarLogN     <- c(0,rep(1,9))
#- Set observation variance parameters
  #lapply(indices[1:2],function(x){dat <- x@index[,drop=T];return(apply(dat,1,sd,na.rm=T)/apply(dat,1,mean,na.rm=T))})
  #apply(catch.n(stock)[,drop=T],1,sd,na.rm=T)/apply(catch.n(stock)[,drop=T],1,mean,na.rm=T))
conf$keyVarObs[1,]  <- c(0,1,2,2,3,3,4,4,4,4)
conf$keyVarObs[2,]  <- c(0,0,1,2,3,3,rep(NA,4))   + max(conf$keyVarObs[1,],na.rm=T) + 1
conf$keyVarObs[3,]  <- c(0,0,0,1,2,3,3,rep(NA,3)) + max(conf$keyVarObs[2,],na.rm=T) + 1
conf$keyVarObs[4,]  <- c(0,rep(NA,9))             + max(conf$keyVarObs[3,],na.rm=T) + 1
if(run == "LPUE_UK_add")
  conf$keyVarObs[5,]<- c(0,rep(NA,9))             + max(conf$keyVarObs[4,],na.rm=T) + 1

#- Set observation correlation structure
conf$obsCorStruct[] <- c("ID","AR","AR","ID")
if(run == "LPUE_UK_add")
  conf$obsCorStruct[] <- c("ID","AR","ID","ID","ID")
conf$keyBiomassTreat[4:5] <- 1

#conf$keyCorObs[2,]  <- c(rep(0,5),rep(NA,4))
conf$keyCorObs[2,]  <- c(0,rep(1,4),rep(NA,4))

#- Set Fbar range
conf$fbarRange      <- range(stock)[c("minfbar","maxfbar")]
conf                <- setNAtominone(conf)

#- Create default parameter settings
par                 <- defpar(dat,conf)
if(run == "trim89" | sens == "trim89"){
  load(file.path(outPath,"pin_trim98_assessmentOut.RData"))
  par$logN[]        <- pin.logN
  par$logF[]        <- pin.logF
  par$missing       <- pin.missing
}


### ------------------------------------------------------------------------------------------------------
###   4. Run assessment
### ------------------------------------------------------------------------------------------------------

fit                 <- sam.fit(dat,conf,par)
pin.logN            <- fit$pl$logN; pin.logF <- fit$pl$logF; pin.missing <- fit$pl$missing
forecastfit         <- try(forecast(fit, fscale=c(1,1,1)))
retro               <- try(retro(fit,7))
LO                  <- try(leaveout(fit))
resids              <- try(residuals(fit))
save(fit,forecastfit,retro,LO,resids,dat,conf,stock,indices,par,file=file.path(outPath,paste0(run,"_",sens,"assessmentOut.RData")))

retro <- NULL
LO    <- NULL

### ------------------------------------------------------------------------------------------------------
###   5. Diagnostics
### ------------------------------------------------------------------------------------------------------
pdf(file.path(outPath,paste0(run,"_",sens,"assessmentOut.pdf")))

for(i in 1:(length(surveys)+1))
  fitplot(fit,fleets=i)
plot(resids)
obscorrplot(fit)
srplot(fit)
ssbplot(fit)
fbarplot(fit)
recplot(fit)
ssbplot(retro)
fbarplot(retro)
recplot(retro)
ssbplot(LO)
fbarplot(LO)
recplot(LO)
ssbplot(forecastfit)
fbarplot(forecastfit)
recplot(forecastfit)

catchdat <- data.frame(value=exp(fit$pl$logFpar),survey=rep(names(surveys),times=unlist(lapply(apply(setminonetoNA(conf)$keyLogFpar[-1,],1,function(x){return(na.omit(unique(x)))}),length))),param=1:length(exp(fit$pl$logFpar)))
for(i in 2:nrow(catchdat))
  catchdat$param[i] <- ifelse(catchdat$survey[i-1] == catchdat$survey[i],catchdat$param[i-1]+1,1)

xyplot(value ~ param|survey,data=catchdat,type="b",scales=list(y="free"),xlab="Parameters",ylab="Catchability")

barplot(exp(fit$pl$logSdLogObs),names=rep(c("Catch",names(surveys)),times=unlist(lapply(apply(setminonetoNA(conf)$keyVarObs,1,function(x){return(na.omit(unique(x)))}),length))),
        col=an(as.factor(rep(c("Catch",names(surveys)),times=unlist(lapply(apply(setminonetoNA(conf)$keyVarObs,1,function(x){return(na.omit(unique(x)))}),length))))),las=2,ylab="Observation variance")

plot(exp(fit$pl$logSdLogFsta))

dev.off()

### ------------------------------------------------------------------------------------------------------
###   6. Further analyses
### ------------------------------------------------------------------------------------------------------

#- Further finetuning
fit3a <- runwithout(fit,year=1996,      fleet=4)
fit3b <- runwithout(fit,year=1996:1997, fleet=rep(4,2))
fit3c <- runwithout(fit,year=1996:1998, fleet=rep(4,3))



