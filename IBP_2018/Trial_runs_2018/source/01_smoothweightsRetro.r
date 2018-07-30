#-------------------------------------------------------------------------------
# IBPNEW 2012
# Code to generate the .dat file for input into the ADMB model
# Uses Lowestoft format input files
# David Miller, some code adapted from FLR-project
# Date: 3 Oct 2012
#-------------------------------------------------------------------------------

## Source methods/functions
source(paste(codePath,"nsea_functions.r",sep=""))

### ------------------------------------------------------------------------------------------------------
###  1. Settings
### ------------------------------------------------------------------------------------------------------


# Units
Units <- c("tonnes","thousands","kg")     # totals (product of the next two), numbers, weights-at-age

# Last Age  (first set at 1)
maxA <- 10
# Is this a plusgroup?
pGrp <- T

# Number of knots in annual F spline (survey/catch splines have 4 knots by default)
Fknots <- 5

### ------------------------------------------------------------------------------------------------------
###   2. Read and process assessment input data
### ------------------------------------------------------------------------------------------------------
## Read stock data
stock              <- readFLStock(paste(dataPath,"index_raw.txt", sep=""))
stock              <- window(stock,end=iYr)
units(stock)[1:17] <- as.list(c(rep(Units,4), "NA", "NA", "f", "NA", "NA"))
stkWtPgrp <- stock.wt(stock)[maxA,]
if (pGrp) stock <- setPlusGroup(stock, plusgroup=maxA) else stock <- trim(stock, age=1:maxA)
stock.wt(stock)[maxA,] <- stkWtPgrp 

# Number of years
years <- as.numeric(range(stock)["minyear"]:range(stock)["maxyear"])
numYr <- as.numeric(range(stock)["maxyear"]-range(stock)["minyear"]+1)
# Number of ages
numAges <- length(1:maxA)

startyr <- range(stock)[["minyear"]]## Read index data

### ------------------------------------------------------------------------------------------------------
###   3. Calculate splines
### ------------------------------------------------------------------------------------------------------
# need to add decimal places if value is an integer (otherwise ADMB reads in the rest as integers too) - so replace 1 or 0 with "1.00" or "0.00"
Fspline <- t(matrix(bs(1:numYr,df=Fknots,intercept=T),ncol=Fknots))
Fspline[Fspline==0|Fspline==1] <- format(Fspline[Fspline==0|Fspline==1],nsmall=2)

### ------------------------------------------------------------------------------------------------------
###   5. Create .dat file
### ------------------------------------------------------------------------------------------------------
turbDAT <- function(stock, numYr, numAges, Fspline){
  cat("#############\n")
  cat("# ",name(stock),"\n")
  cat("# Created:",format(Sys.time(), "%d%b%Y_%Hh%M"),"\n")
  cat("# years:",range(stock)["minyear"],"-",range(stock)["maxyear"]," ; ages:",range(stock)["min"],"-",range(stock)["max"], "; number of knots in time spline \n")
  cat(numYr,numAges,Fknots,"\n")
  cat("#############\n")
  quants <- mcf(list(round(landings.wt(stock),3),round(stock.wt(stock),3)))
  tquants <- lapply(quants,function(x){x <- matrix(x,nrow=dim(x)[1]);t(x);})
  for (ii in 1:length(tquants)){
    tquants[[ii]][is.na(tquants[[ii]])] <- round(-1,0)
    write.table(tquants[[ii]], row.names=F, col.names=F,quote=F)
    cat("\n");
  }
  cat("#############\n")  
  cat("# Annual F spline","\n")
  write.table(Fspline, row.names=F, col.names=F,quote=F)

} # end of function

capture.output(turbDAT(stock, numYr, numAges,  Fspline), file=paste(codePath,"turWsmooth\\Turbot_2014", if(pGrp) "_PG.dat" else "_noPG.dat",sep=""))

# run assessment and print cwt and swt
setwd(paste(codePath,"turWsmooth\\", sep=""))
system("turbot_2014_PG.exe")
report <- readLines("turbot_2014_PG.rep")

#look at smoothed weights
(report[which(report=="Estimated CWT") +(1:numYr)])
(report[which(report=="Estimated SWT") +(1:numYr)])

#plot observed and smoothed stock weights and catch weights
sswts <- matrix(as.numeric(unlist(strsplit(report[which(report=="Estimated SWT") +(1:numYr)], " "))), ncol=numYr)[-1,]
scwts <- matrix(as.numeric(unlist(strsplit(report[which(report=="Estimated CWT") +(1:numYr)], " "))), ncol=numYr)[-1,]

stock@stock.wt[,ac(1981:iYr)] <- sswts[,-c(1:6)]
stock@catch.wt[,ac(1981:iYr)] <- scwts[,-c(1:6)]
stock@landings.wt[,ac(1981:iYr)] <- scwts[,-c(1:6)]


sswts2017 <- sswts
scwts2017 <- scwts

sswts2016 <- sswts
scwts2016 <- t(read.table("D:/Repository/Turbot/assessment runs/lowestoft files/weca.txt",skip=5))

sswts2015 <- sswts
scwts2015 <- scwts

sswts2014 <- sswts
scwts2014 <- scwts

sswts2013 <- sswts
scwts2013 <- scwts

sswts2012 <- sswts
scwts2012 <- scwts

sswts2011 <- sswts
scwts2011 <- scwts

sswts2010 <- sswts
scwts2010 <- scwts

save(sswts2017,sswts2016,sswts2015,sswts2014,sswts2013,sswts2012,sswts2011,sswts2010,file=file.path(outPath,"sswtsRetro.RData"))
save(scwts2017,scwts2016,scwts2015,scwts2014,scwts2013,scwts2012,scwts2011,scwts2010,file=file.path(outPath,"scwtsRetro.RData"))

par(mfrow=c(3,4))
for(iAge in 1:10){
  for(iYr in 2017:2010){
    assign("dat",get(paste0("sswts",iYr)))
    yrange <- range(dat[iAge,]) * c(0.9,1.1)
    if(iYr == 2017)
      plot(dat[iAge,],x=1975:iYr,xlab="Years",ylab="Weight (kg)",main=paste("Stock weigth:",iAge),ylim=yrange,las=1,pch=substr(iYr,4,4),type="b",col=abs(an(substr(iYr,4,4))+1),cex=0.75)
    if(iYr != 2017)
      points(dat[iAge,],x=1975:iYr,pch=substr(iYr,4,4),type="b",col=an(substr(iYr,4,4)),cex=0.75)
  }
}