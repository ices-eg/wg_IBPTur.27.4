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



### ------------------------------------------------------------------------------------------------------
###   2. make new plots
### ------------------------------------------------------------------------------------------------------

#Landings
plot(x=startyr:(assYear-1), y=landings(stock)/1000, xlim=c(1957,(assYear-1)), ylim=c(0,1.1*max(landings(stock)))/1000, type="l", xlab="Year",lty=1, ylab= "Landings ('000 t)",main="Landings", las=1, lwd=2, col="blue")
grid()

#Landings and discards (manually added discards, which are in a diton file in lowestoft files and come from intercatch)
plot(x=startyr:(assYear-1), y=landings(stock)/1000, xlim=c(1957,(assYear-1)), ylim=c(0,1.1*max(landings(stock)))/1000, type="l", xlab="Year",lty=1, ylab= "Landings ('000 t)",main="Landings", las=1, lwd=2, col="blue")
lines(x=2013:2016, y=c(97,158,112,666)/1000, col="red", lwd=2)
grid()


bsize <- 0.08 
plc <-  (landings.n(stock))[,ac(1975: (assYear-1))]
# bubbles(age~year, data=resL, col=c("black","black"), bub.scale=10, pch=c(21,21), fill=resL>0, xlim=c(1957,(assYear-1)), ylim=c(0,(maxA+1)), ylab= "Age", main="Landings residuals")
nmages <-  length(dimnames(plc)[[1]])
ylims <- 0.5*nmages+0.25
plot(NA,NA,main="Landings at age", xlab="Year", ylab="Age",xlim=c(as.numeric(min(dimnames(plc)[[2]])),as.numeric(max(dimnames(plc)[[2]]))), yaxt="n", ylim=c(0,ylims))
axis(2,at=seq(0.5,0.5*nmages,by=0.5),labels=as.numeric(dimnames(plc)[[1]]))
for (i in as.numeric(dimnames(plc)[[1]])) {
  radius <- as.numeric(sqrt(abs(plc[i,]) / pi)) *bsize
  points(dimnames(plc)[[2]], rep(0.5*i,ncol(plc)),cex=radius*2, pch=21, col=c("black", "blue")[1+as.numeric(plc[i,]>0)], bg=alpha(c("black", "blue")[1+as.numeric(plc[i,]>0)],0.5))
}
text(2005,0.1,paste("min = ",round(min(c(plc)[!is.infinite(c(plc))],na.rm=T),2),"; max = ",round(max(c(plc)[!is.infinite(c(plc))],na.rm=T),2) ,sep=""), cex=1, pos=4)

#z from catch curves (based on landings)

#Surveys
plot(indices[[1]])
plot(indices[[2]])
plot(indices[[3]])


#### Run indices diagnostics
idxcrop <-  indices
indsN01 <- FLQuants(lapply( mcf( lapply(idxcrop, index)), function(x){x <- FLQuant(aperm(apply(x@.Data, c(1,3,4,5,6), scale),c(2,1,3,4,5,6)), dimnames= dimnames(x))}))

names(indsN01)   <- names(indices)
akey             <- simpleKey(text=names(indsN01), points=F, lines=F, columns=2, cex=1.5, col=c("red","black","blue","gray","orange","magenta"))
#akey$text$lab[1] <- "BTS-ISIS"; akey$text$lab[2] <- "SNS"
xyplot(data~year | factor(age), data=indsN01, type="b", key=akey, groups=qname, pch=19, 
       col=c("red","black","blue","gray","orange","magenta"),as.table=TRUE, scales="free",  layout=c(3,3), xlim=c(1995,(assYear-1)), ylim=c(-1.5,3.5))

for (i in c(1,2,3)) {
  plotInternalConsistency(indices[[i]], mark.last=T)
}

#plot full and survey indices and those that are cut from when ageing started
#SNS
plotInternalConsistency(indices[[1]], mark.last=T)
plotInternalConsistency(window(indices[[1]], start=2004), mark.last=T)

#BTS
plotInternalConsistency(indices[[2]], mark.last=T)
plotInternalConsistency(window(indices[[2]],start=1996), mark.last=T)


#plot unstructured indices
indices_unstruct  <- readFLIndices(paste(dataPath,"fleet_tmb_longer_LPUE.txt", sep=""), na.strings="-1")

plot(index(indices_unstruct[[3]]), ylim=c(0,110))


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

