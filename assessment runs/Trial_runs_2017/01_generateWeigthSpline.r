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
library(FLAssess); library(stockassessment)
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
stock_Name      <- "tur-nsea"
# Year (= year when assessment conducted.  i.e. have data up to assYear-1)
assYear         <- 2017
retroYr         <- 2017
endYear         <- min(assYear,retroYr)-1
Units           <- c("tonnes","thousands","kg")     # totals (product of the next two), numbers, weights-at-age
maxA            <- 10
pGrp            <- T
minFbar         <- 2
maxFbar         <- 6

### ------------------------------------------------------------------------------------------------------
###   2. Read and process assessment input data
### ------------------------------------------------------------------------------------------------------

## Read stock data
stock               <- readFLStock(paste(dataPath,"index_raw.txt", sep=""))
units(stock)[1:17]  <- as.list(c(rep(Units,4), "NA", "NA", "f", "NA", "NA"))
# Change old ages landings to 0 (instead of NA) to allow plusgroup weight calculation
# Calculate plusgroup weight for stock
stkWtPgrp           <- stock.wt(stock)[maxA,]
# Apply max age
if (pGrp) stock     <- setPlusGroup(stock, plusgroup=maxA) else stock <- trim(stock, age=1:maxA)
# Trim if doing retro
stock.wt(stock)[maxA,] <- stkWtPgrp 

# Number of years
years               <- as.numeric(range(stock)["minyear"]:range(stock)["maxyear"])
numYr               <- as.numeric(range(stock)["maxyear"]-range(stock)["minyear"]+1)
# Number of ages
numAges             <- length(1:maxA)
startyr             <- range(stock)[["minyear"]]## Read index data

# NOTE: all index values should run to (assYear-1) - because only #yrs read in, so start of index is calculated back from (assYear-1) and this.
indices             <- readFLIndices(paste(dataPath,"fleet_tmb_longer_LPUE.txt", sep=""), na.strings="-1")
indices             <- FLIndices(list(indices[[1]],indices[[2]],indices[[4]]))
indexVals           <- lapply(indices, index)
numIndices          <- length(indexVals)
indMPs              <- list()
for(ind in names(indices))
  indMPs[[ind]]     <- as.numeric((range(indices[[ind]])["startf"]+range(indices[[ind]])["endf"])/2)


#### Plot observed stock weights
maxA  <- 10
cols <- rich.colors(maxA)
stock.wt(stock)[stock.wt(stock)==0] <- NA
plot(x=1975:(assYear-1), y=stock.wt(stock)[1,], ylim=c(0,9.1), xlab="Year", cex.lab=1.2,ylab="Weights (kg)", main="Stock weight@age", lwd=2, col="white" )
for (aa in 1:maxA){
  points(x=1975:(assYear-1),y=as.numeric(stock.wt(stock)[aa,]), pch= if (aa<10) as.character(aa) else "+", col=cols[aa])
}

landings.wt(stock)[landings.wt(stock)==0] <- NA
##----landing---
plot(x=1975:(assYear-1), y=landings.wt(stock)[1,], ylim=c(0,9.1), xlab="Year", cex.lab=1.2,ylab="Weights (kg)", main="Landing weight@age", lwd=2, col="white" )
for (aa in 1:maxA){
  points(x=1975:(assYear-1),y=as.numeric(landings.wt(stock)[aa,]), pch= if (aa<10) as.character(aa) else "+", col=cols[aa])
}


discards.wt(stock)[discards.wt(stock) < 0.05] <- NA
##----discards---
plot(x=1975:(assYear-1), y=discards.wt(stock)[1,], ylim=c(0,1.1), xlab="Year", cex.lab=1.2,ylab="Weights (kg)", main="Discards weight@age", lwd=2, col="white" )
for (aa in 1:8){
  points(x=1975:(assYear-1),y=as.numeric(discards.wt(stock)[aa,]), pch= if (aa<10) as.character(aa) else "+", col=cols[aa])
}

catch.wt(stock)[discards.wt(stock) < 0.05] <- NA
##----catch---
plot(x=1975:(assYear-1), y=catch.wt(stock)[1,], ylim=c(0,9.1), xlab="Year", cex.lab=1.2,ylab="Weights (kg)", main="Catch weight@age", lwd=2, col="white" )
for (aa in 1:10){
  points(x=1975:(assYear-1),y=as.numeric(catch.wt(stock)[aa,]), pch= if (aa<10) as.character(aa) else "+", col=cols[aa])
}

### ------------------------------------------------------------------------------------------------------
###   4. Calculate splines
### ------------------------------------------------------------------------------------------------------
# Number of knots in annual F spline (survey/catch splines have 4 knots by default)
Fknots            <- 5
numFplateau       <- 1  #if 1: out to age 7, if 2: one to age 7 and one to age 6  (in that order)
# need to add decimal places if value is an integer (otherwise ADMB reads in the rest as integers too) - so replace 1 or 0 with "1.00" or "0.00"
selSplines        <- list()
# survey selectivity
selSplines[[1]]   <- t(matrix(bs(1:7,4,T),ncol=4))
selSplines[[1]][selSplines[[1]]==0|selSplines[[1]]==1] <- format(selSplines[[1]][selSplines[[1]]==0|selSplines[[1]]==1],nsmall=2)
# survey selectivity
selSplines[[2]]   <- t(matrix(bs(1:6,4,T),ncol=4))
selSplines[[2]][selSplines[[2]]==0|selSplines[[2]]==1] <- format(selSplines[[2]][selSplines[[2]]==0|selSplines[[2]]==1],nsmall=2)
# Annual F
Fspline           <- t(matrix(bs(1:numYr,df=Fknots,intercept=T),ncol=Fknots))
Fspline[Fspline==0|Fspline==1] <- format(Fspline[Fspline==0|Fspline==1],nsmall=2)

### ------------------------------------------------------------------------------------------------------
###   5. Create .dat file
### ------------------------------------------------------------------------------------------------------

turbDAT <- function(stock, numYr, numAges, indexVals, indMPs, selSplines, Fspline){
cat("#############\n")
cat("# ",name(stock),"\n")
cat("# Created:",format(Sys.time(), "%d%b%Y_%Hh%M"),"\n")
cat("# years:",range(stock)["minyear"],"-",range(stock)["maxyear"]," ; ages:",range(stock)["min"],"-",range(stock)["max"],"; Fbar range", minFbar,"-", maxFbar, "; number of knots in time spline \n")
cat(numYr,numAges, " 7 ", minFbar, maxFbar, " 4 ",Fknots,"\n")
cat(indMPs[["SNS"]],indMPs[["BTS-ISIS"]],indMPs[["Dutch_BT2_LPUE_ages"]],"\n")
cat("#############\n")
    
quants <- mcf(list(landings.n(stock),round(landings.wt(stock),3),round(stock.wt(stock),3), indexVals[["SNS"]],indexVals[["BTS-ISIS"]],indexVals[["Dutch_BT2_LPUE_ages"]]))
tquants <- lapply(quants,function(x){x <- matrix(x,nrow=dim(x)[1]);t(x);})
#lapply(tquants,function(x){x[is.na(x)] <- round(-1,0) ; write.table(x, row.names=F, col.names=F,quote=F);cat("\n");})

for (ii in 1:length(tquants)){
  if (!(ii %in% c(2,3))){ tquants[[ii]] <- tquants[[ii]] + min(tquants[[ii]][!tquants[[ii]]==0], na.rm=T)/2}
  tquants[[ii]][is.na(tquants[[ii]])] <- round(-1,0)
  write.table(tquants[[ii]], row.names=F, col.names=F,quote=F)
  cat("\n");
}

cat("#############\n")
cat("# landings","\n")
cat(landings(stock),"\n")
cat("#############\n")
cat("# M","\n")
cat(m(stock)[1,1],"\n")
cat("#############\n")
cat("# Maturity (females)","\n")
tmp <- matrix(mat(stock)[,1]@.Data,ncol=1)
tmp[is.na(tmp)] <- round(-1,0)
write.table(tmp, row.names=F, col.names=F,quote=F)
#cat("\n")
cat("#############\n")
cat("# Selectivity spline (surveys): 4 knot, 7 ages (last ages are equal)","\n")
write.table(selSplines[[1]], row.names=F, col.names=F,quote=F)
#cat("\n")
if (numFplateau>1) {
cat("#############\n")
cat("# Selectivity spline (catch, LPUE): 4 knot, 6 ages (last ages are equal)","\n")
write.table(selSplines[[2]], row.names=F, col.names=F,quote=F)
}
##cat("\n")
cat("#############\n")  
cat("# Annual F spline","\n")
write.table(Fspline, row.names=F, col.names=F,quote=F)
#cat(,"\n")

} # end of function

capture.output(turbDAT(stock, numYr, numAges, indexVals, indMPs, selSplines, Fspline), file=paste(outPath,"Exploratory assessment VBLG PG\\Turbot_2014", if(pGrp) "_PG.dat" else "_noPG.dat",sep=""))

# run assessment and print cwt and swt
setwd(paste(outPath,"Exploratory assessment VBLG PG\\", sep=""))
system("turbot_2014_PG.exe")
report <- readLines("turbot_2014_PG.rep")

#look at smoothed weights
(report[which(report=="Estimated CWT") +(1:numYr)])
(report[which(report=="Estimated SWT") +(1:numYr)])


#plot observed and smoothed stock weights and catch weights
sswts <- matrix(as.numeric(unlist(strsplit(report[which(report=="Estimated SWT") +(1:numYr)], " "))), ncol=numYr)[-1,]
scwts <- matrix(as.numeric(unlist(strsplit(report[which(report=="Estimated CWT") +(1:numYr)], " "))), ncol=numYr)[-1,]

save(sswts,scwts,file=file.path(outPath,"smoothedWeights.RData"))

maxA  <- 10
cols <- rich.colors(maxA)
stock.wt(stock)[stock.wt(stock)==0] <- NA
plot(x=1975:(assYear-1), y=stock.wt(stock)[1,], ylim=c(0,10.1), xlab="Year", cex.lab=1.2,ylab="Weights (kg)", main="Stock weight@age", lwd=2, col="white" )
for (aa in 1:maxA){
  points(x=1975:(assYear-1),y=as.numeric(stock.wt(stock)[aa,]), pch= if (aa<10) as.character(aa) else "+", col=cols[aa])
  lines(x=1975:(assYear-1), y=sswts[aa,], col=cols[aa])
  }

maxA  <- 10
cols <- rich.colors(maxA)
landings.wt(stock)[landings.wt(stock)==0] <- NA
plot(x=1975:(assYear-1), y=landings.wt(stock)[1,], ylim=c(0,10.1), xlab="Year", cex.lab=1.2,ylab="Weights (kg)", main="Landings weight@age", lwd=2, col="white" )
for (aa in 1:maxA){
  points(x=1975:(assYear-1),y=as.numeric(landings.wt(stock)[aa,]), pch= if (aa<10) as.character(aa) else "+", col=cols[aa])
  lines(x=1975:(assYear-1), y=scwts[aa,], col=cols[aa])
}


btsagesamp <- read.sas7bdat( "w://IMARES/Data/ICES-WG/Demersale werkgroep WGNSSK/2017/Stock/tur-nsea/turbot_survey_cpue/tur_bri bts/age_tur.sas7bdat")
head(btsagesamp)
btsagesamp$age <- btsagesamp$year - btsagesamp$yc  
btsagesamp <- btsagesamp[!(is.na(btsagesamp$weight) | is.na(btsagesamp$yc)),]


gambtswts <- gam(weight ~ te(yc,age),data=btsagesamp)
summary(gambtswts)
length(gambtswts$fitted)
nrow(btsagesamp)
btsagesamp$fitted <- gambtswts$fitted


aggbtswts <- aggregate(weight ~ yc + age + year, data=btsagesamp, FUN="mean")

par(mfrow=c(1,1))
plot((weight/1000)~year, data=aggbtswts[aggbtswts$age==1,],xlim=c(1975,assYear-1),ylim=c(0,9.1), xlab="Year", cex.lab=1.2,ylab="Weights (kg)", main="BTS weight@age", lwd=2, col="white" )
for (aa in 1:maxA){
  points((weight/1000)~year, data=aggbtswts[aggbtswts$age==aa,],pch= if (aa<10) as.character(aa) else "+", lwd=2, col=cols[aa])
#  lines(x=1975:(assYear-1), y=scwts[aa,], col=cols[aa])
}

par(mfrow=c(2,5), mar=c(4.1,4.1,0.2,0.2))
for (aa in 1:maxA){
  plot((weight/1000)~year, data=aggbtswts[aggbtswts$age==aa,],pch= 15, lwd=2, col=cols[aa], ylim=c(0,9.1), xlim=c(1975,2017))
  #  lines(x=1975:(assYear-1), y=scwts[aa,], col=cols[aa])
  points(x=1975:(assYear-1),y=as.numeric(stock.wt(stock)[aa,]), pch= 17, col=cols[aa])
}

##look at sop

scwts <- matrix(as.numeric(unlist(strsplit(report[which(report=="Estimated CWT") +(1:numYr)], " "))), ncol=numYr)[-1,]

#SOP?
plot.default(x=years, y=landings(stock),ylim=c(0,8000), type="l" )
temp <- apply(landings.n(stock)[,,drop=T]*scwts,2,sum,na.rm=T)
temp[temp<1] <- NA
lines(x=years, y=temp, lty=2)

sop <- temp/landings(stock)

landings.n(stock) <- sweep(landings.n(stock),2,sop,"/")

