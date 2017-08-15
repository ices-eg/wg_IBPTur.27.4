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
## Stock name
stock_Name      <- "tur-nsea"

## Assessment settings
# Year (= year when assessment conducted.  i.e. have data up to assYear-1)
assYear <- 2017
retroYr <- 2017
# last year of data:
endYear <- min(assYear,retroYr)-1  

# Units
Units <- c("tonnes","thousands","kg")     # totals (product of the next two), numbers, weights-at-age

# Last Age  (first set at 1)
maxA <- 10
# Is this a plusgroup?
pGrp <- T

# number of F plateaus
numFplateau <- 1  #if 1: out to age 7, if 2: one to age 7 and one to age 6  (in that order)

# Number of knots in annual F spline (survey/catch splines have 4 knots by default)
Fknots <- 5

#Fbar age range
minFbar <- 2
maxFbar <- 6
### ------------------------------------------------------------------------------------------------------
###   2. Read and process assessment input data
### ------------------------------------------------------------------------------------------------------
## Read stock data
stock              <- readFLStock(paste(dataPath,"index_raw.txt", sep=""))
units(stock)[1:17] <- as.list(c(rep(Units,4), "NA", "NA", "f", "NA", "NA"))
# Change old ages landings to 0 (instead of NA) to allow plusgroup weight calculation
# Calculate plusgroup weight for stock
stkWtPgrp <- stock.wt(stock)[maxA,]
# Apply max age
if (pGrp) stock <- setPlusGroup(stock, plusgroup=maxA) else stock <- trim(stock, age=1:maxA)
# Trim if doing retro
stock.wt(stock)[maxA,] <- stkWtPgrp 

# Number of years
years <- as.numeric(range(stock)["minyear"]:range(stock)["maxyear"])
numYr <- as.numeric(range(stock)["maxyear"]-range(stock)["minyear"]+1)
# Number of ages
numAges <- length(1:maxA)

startyr <- range(stock)[["minyear"]]## Read index data


# NOTE: all index values should run to (assYear-1) - because only #yrs read in, so start of index is calculated back from (assYear-1) and this.
indices            <- readFLIndices(paste(dataPath,"fleet_tmb_longer_LPUE.txt", sep=""), na.strings="-1")

indices <- FLIndices(list(indices[[1]],indices[[2]],indices[[4]]))

indexVals <- lapply(indices, index)
numIndices <- length(indexVals)
indMPs <- list()
for (ind in names(indices)) indMPs[[ind]] <- as.numeric((range(indices[[ind]])["startf"]+range(indices[[ind]])["endf"])/2)



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

### ------------------------------------------------------------------------------------------------------
###   3. Calculate splines
### ------------------------------------------------------------------------------------------------------
# need to add decimal places if value is an integer (otherwise ADMB reads in the rest as integers too) - so replace 1 or 0 with "1.00" or "0.00"
selSplines <- list()
# survey selectivity
selSplines[[1]] <- t(matrix(bs(1:7,4,T),ncol=4))
selSplines[[1]][selSplines[[1]]==0|selSplines[[1]]==1] <- format(selSplines[[1]][selSplines[[1]]==0|selSplines[[1]]==1],nsmall=2)
# survey selectivity
selSplines[[2]] <- t(matrix(bs(1:6,4,T),ncol=4))
selSplines[[2]][selSplines[[2]]==0|selSplines[[2]]==1] <- format(selSplines[[2]][selSplines[[2]]==0|selSplines[[2]]==1],nsmall=2)
# Annual F
Fspline <- t(matrix(bs(1:numYr,df=Fknots,intercept=T),ncol=Fknots))
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

#try effect of SOP 



setwd("d:\\tur-nsea\\lowestoft files\\")
cn<-read.ices("canum.txt")
cn <- t(landings.n(stock)[,,drop=T])
cw<-read.ices("weca.txt")
dw<-read.ices("weca.txt")
lw<-read.ices("weca.txt")
mo<-read.ices("matprop.txt")
nm<-read.ices("natmor.txt")
pf<-read.ices("fprop.txt")
pm<-read.ices("mprop.txt")
sw<-read.ices("west.txt")
lf<-read.ices("lf.txt")

surveys<-read.ices("fleet_tmb_longer_LPUE.txt")

dat<-setup.sam.data(surveys=surveys,
                    residual.fleet=cn, 
                    prop.mature=mo, 
                    stock.mean.weight=sw, 
                    catch.mean.weight=cw, 
                    dis.mean.weight=dw, 
                    land.mean.weight=lw,
                    prop.f=pf, 
                    prop.m=pm, 
                    natural.mortality=nm, 
                    land.frac=lf)

conf<-defcon(dat)
conf$keyLogFsta[1,] <- c(0,    1,    2,    3,    4,    5,    6,    7,    8,     8)

conf$keyLogFpar[2,] <- c(0,    1 ,   2  ,  3 ,   4 ,   4  ,  4,   -1,   -1,    -1)
conf$keyLogFpar[3,] <- c(5,    6 ,   7  ,  8 ,   9 ,   9  ,  9,   -1,   -1,    -1)
conf$keyLogFpar[4,] <- c(10,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,    -1)

conf$keyVarF[1,] <- c(0, 1, 2, 2, 2, 2, 2, 2, 2, 2)

conf$keyVarLogN <-  c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1)

conf$keyVarObs[1,]<-c(0,1,1,1,1,1,1,1,1,1)
conf$keyVarObs[2,]<-c(2,2,2,2,2,2,2,-1,-1,-1)
conf$keyVarObs[3,]<-c(3,3,3,3,3,3,3,-1,-1,-1)
conf$keyVarObs[4,]<-c(4,-1,-1,-1,-1,-1,-1,-1,-1,-1)

conf$obsCorStruct[]<-c("ID","AR","AR","ID")
conf$keyCorObs[2,]<-c(0,0,0,0,0,0,-1,-1,-1)
conf$keyCorObs[3,]<-c(1,1,1,1,1,1,-1,-1,-1)
conf$fbarRange<-c(2,6)
conf$corFlag<-2
par<-defpar(dat,conf)
fit<-sam.fit(dat,conf,par)
ssbplot(fit)
AIC(fit)
ssbplot(fit, xlim=c(1975,2018))
points(x=2017,y=forecast(fit,fscale=c(1,1 ))[2,7])

R <- summary(fit)[,1]
SSB <- summary(fit)[,4]
plot(x=SSB[-length(SSB)],y=R[-1], type="l", xlim=c(0,15000), ylim=c(0,7000))
abline(lm(R[-1]~ SSB[-length(SSB)]))
srplot(fit)
fitme <- glm(log(R[-1])~ SSB[-length(SSB)])
lines(x=SSB[-length(SSB)], y=exp(fitme$fitted.values), col="blue", lty=2)
abline( glm(R[-1]~ SSB[-length(SSB)]), col="blue", lty=2)


fit$sdrep
fval <- exp(fit$pl$logF)
fval <- rbind(fval,fval[9,])

Ns <- exp(fit$pl$logN)

survivors <- Ns*exp(-(fval + m(stock)[,,drop=T]@.Data))
survivors[9,] <- survivors[9,]  + survivors[10,] 
survivors <- survivors[,ncol(survivors)]

LO <- leaveout(fit)
ssbplot(LO)


surveys<-read.ices("fleet_tmb_longer_LPUE_SNS1-3.txt")

dat<-setup.sam.data(surveys=surveys,
                    residual.fleet=cn, 
                    prop.mature=mo, 
                    stock.mean.weight=sw, 
                    catch.mean.weight=cw, 
                    dis.mean.weight=dw, 
                    land.mean.weight=lw,
                    prop.f=pf, 
                    prop.m=pm, 
                    natural.mortality=nm, 
                    land.frac=lf)

conf<-defcon(dat)
conf$keyLogFsta[1,] <- c(0,    1,    2,    3,    4,    5,    6,    7,    7,     7)

conf$keyLogFpar[2,] <- c(0,    1 ,   2  ,  -1 ,   -1 ,   -1  ,  -1,   -1,   -1,    -1)
conf$keyLogFpar[3,] <- c(3,    4 ,   5  ,  6 ,   7 ,   7  ,  7,   -1,   -1,    -1)
conf$keyLogFpar[4,] <- c(8,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,    -1)

conf$keyVarF[1,] <- c(0, 1, 2, 2, 2, 2, 2, 2, 2, 2)

conf$keyVarLogN <-  c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1)

conf$keyVarObs[1,]<-c(0,1,1,1,1,1,1,1,1,1)
conf$keyVarObs[2,]<-c(2,2,2,-1,-1,-1,-1,-1,-1,-1)
conf$keyVarObs[3,]<-c(3,3,3,3,3,3,3,-1,-1,-1)
conf$keyVarObs[4,]<-c(4,-1,-1,-1,-1,-1,-1,-1,-1,-1)

conf$obsCorStruct[]<-c("ID","AR","AR","ID")
conf$keyCorObs[2,]<-c(0,0,0,-1,-1,-1,-1,-1,-1)
conf$keyCorObs[3,]<-c(1,1,1,1,1,1,-1,-1,-1)
conf$fbarRange<-c(2,6)
conf$corFlag<-2
par<-defpar(dat,conf)
fit<-sam.fit(dat,conf,par)
AIC(fit)
ssbplot(fit)
fbarplot(fit)
forecastsam <- forecast(fit, fscale=c(1,1))
ssbplot(forecastsam)

#allow slightly more flexibility in F
conf$keyLogFsta[1,] <- c(0,    1,    2,    3,    4,    5,    6,    7,    8,     8)
par<-defpar(dat,conf)
fit2<-sam.fit(dat,conf,par)
AIC(fit2)
forecastsam <- forecast(fit2, fscale=c(1,1))
ssbplot(forecastsam)

res=residuals(fit2)
plot(res)
obscorrplot(fit2)
fitplot(fit2,fleets=1, ylim=c(-2,8))
fitplot(fit2,fleets=2, ylim=c(-3,7))
fitplot(fit2,fleets=3, ylim=c(-8,2))


#true zeroes for surveys
surveys<-read.ices("fleet_tmb_longer_LPUE_SNS1-3_true_zeroes.txt")
dat<-setup.sam.data(surveys=surveys,
                    residual.fleet=cn, 
                    prop.mature=mo, 
                    stock.mean.weight=sw, 
                    catch.mean.weight=cw, 
                    dis.mean.weight=dw, 
                    land.mean.weight=lw,
                    prop.f=pf, 
                    prop.m=pm, 
                    natural.mortality=nm, 
                    land.frac=lf)

conf<-defcon(dat)
conf$keyLogFsta[1,] <- c(0,    1,    2,    3,    4,    5,    6,    7,    7,     7)

conf$keyLogFpar[2,] <- c(0,    1 ,   2  ,  -1 ,   -1 ,   -1  ,  -1,   -1,   -1,    -1)
conf$keyLogFpar[3,] <- c(3,    4 ,   5  ,  6 ,   7 ,   7  ,  7,   -1,   -1,    -1)
conf$keyLogFpar[4,] <- c(8,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,    -1)

conf$keyVarF[1,] <- c(0, 1, 2, 2, 2, 2, 2, 2, 2, 2)

conf$keyVarLogN <-  c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1)

conf$keyVarObs[1,]<-c(0,1,1,1,1,1,1,1,1,1)
conf$keyVarObs[2,]<-c(2,2,2,-1,-1,-1,-1,-1,-1,-1)
conf$keyVarObs[3,]<-c(3,3,3,3,3,3,3,-1,-1,-1)
conf$keyVarObs[4,]<-c(4,-1,-1,-1,-1,-1,-1,-1,-1,-1)

conf$obsCorStruct[]<-c("ID","AR","AR","ID")
conf$keyCorObs[2,]<-c(0,0,0,-1,-1,-1,-1,-1,-1)
conf$keyCorObs[3,]<-c(1,1,1,1,1,1,-1,-1,-1)
conf$fbarRange<-c(2,6)
conf$corFlag<-2
par<-defpar(dat,conf)

fit3<-sam.fit(dat,conf,par)

AIC(fit3)
ssbplot(fit3)
forecastsam <- forecast(fit3, fscale=c(1,1))
ssbplot(forecastsam)
fbarplot(forecastsam)
res=residuals(fit3)
LO <- leaveout(fit3)
ssbplot(LO)

retro3 <- retro(fit3, year=3)
ssbplot(retro3)


fit3a <- runwithout(fit3,year=1996, fleet=4)
fit3b <- runwithout(fit3,year=1996:1997, fleet=rep(4,2))
fit3c <- runwithout(fit3,year=1996:1998, fleet=rep(4,3))
fit3d <- runwithout(fit3,year=1996:1999, fleet=rep(4,4))
fit3e <- runwithout(fit3,year=1996:2000, fleet=rep(4,5))
fit3f <- runwithout(fit3,year=1996:2001, fleet=rep(4,6))

ssbplot(c(fit3,fit3a, fit3b, fit3c, fit3d, fit3e, fit3f))
ssbplot(fit3f)
AIC(fit3f)
res3f = residuals(fit3f)
plot(res3f)

ssbplot(forecast(fit3f,fscale=c(1,1)))
forecas

#compare previous results to what happens when extending keylogFsta
conf$keyLogFsta[1,] <- c(0,    1,    2,    3,    4,    5,    6,    7,    8,     8)

par<-defpar(dat,conf)

fit4<-sam.fit(dat,conf,par)
fit4a <- runwithout(fit4,year=1996, fleet=4)
fit4b <- runwithout(fit4,year=1996:1997, fleet=rep(4,2))
fit4c <- runwithout(fit4,year=1996:1998, fleet=rep(4,3))
fit4d <- runwithout(fit4,year=1996:1999, fleet=rep(4,4))
fit4e <- runwithout(fit4,year=1996:2000, fleet=rep(4,5))
fit4f <- runwithout(fit4,year=1996:2001, fleet=rep(4,6))

ssbplot(c(fit4,fit4a, fit4b, fit4c, fit4d, fit4e, fit4f))

AIC(fit4)
ssbplot(fit4)
res4= residuals(fit4)
plot(res4)

ssbplot(fit4f)
forecastsam <- forecast(fit4f, fscale=c(1,1))
ssbplot(forecastsam)
