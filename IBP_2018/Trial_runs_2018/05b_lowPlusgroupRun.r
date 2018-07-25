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
outPath   <- paste(Path,"Trial_runs_2018/output/base_pgroup/",sep="")
codePath  <- paste(Path,"Trial_runs_2018/source/",sep="")

## Source methods/functions
source(paste(codePath,"01_smoothweights.r",sep=""))
source(paste(codePath,"03a_setupStockIndices.r",sep=""))

run       <- "base_pgroup"
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
stock               <- TUR
TUR.tun             <- indices
TUR.ctrl            <- FLSAM.control(TUR,TUR.tun)

### ------------------------------------------------------------------------------------------------------
###   3. Setup data structure for SAM assessment
### ------------------------------------------------------------------------------------------------------

TUR.sams  <- list()
TUR.sams.retro <- list()
for(pg in 10:6){
  TUR                 <- setPlusGroup(stock,pg)
  TUR@name            <- paste(TUR@name,"pg",pg,sep=" ")
  TUR.tun             <- indices
  if(pg == 6){
    TUR.tun[["BTS-ISIS"]] <- trim(TUR.tun[["BTS-ISIS"]],age=1:6)
  }
  TUR.ctrl            <- FLSAM.control(TUR,TUR.tun)

  if(pg >9 )
    TUR.ctrl@states["catch unique",]            <- c(0:(pg-3),rep(pg-3,2))
  if(pg <=9)
    TUR.ctrl@states["catch unique",]            <- c(0:(pg-2),(pg-2))
  TUR.ctrl@cor.F                                <- 2
  TUR.ctrl@catchabilities["SNS",ac(1:6)]        <- c(0:2,rep(3,3))         + 101
  if(pg > 6)
    TUR.ctrl@catchabilities["BTS-ISIS",ac(1:7)] <- c(0,0,1,1,rep(2,3))     + 201
  if(pg == 6)
    TUR.ctrl@catchabilities["BTS-ISIS",ac(1:6)] <- c(0,0,1,1,rep(2,2))     + 201
  TUR.ctrl@catchabilities["NL_LPUE",ac(1)]      <- 0                       + 301
#  TUR.ctrl@catchabilities["IBTS_Q1",ac(1)]     <- 0                       + 401
  TUR.ctrl@f.vars["catch unique",]              <- c(0,1,2,2,3,3,3,4,4,4)[1:pg]
  TUR.ctrl@f.vars["catch unique",ncol(TUR.ctrl@f.vars)] <-   TUR.ctrl@f.vars["catch unique",ncol(TUR.ctrl@f.vars)-1]
  TUR.ctrl@logN.vars[]                          <- c(0,rep(1,9))[1:pg]

  TUR.ctrl@obs.vars["catch unique",]            <- c(0,0,1,1,2,2,3,3,4,4)[1:pg]   + 101
  TUR.ctrl@obs.vars["SNS",ac(1:6)]              <- c(0,0,1,1,2,2)           + 201
  if(pg > 6)
    TUR.ctrl@obs.vars["BTS-ISIS",ac(1:7)]       <- c(0,0,1,1,2,2,2)         + 301
  if(pg == 6)
    TUR.ctrl@obs.vars["BTS-ISIS",ac(1:6)]       <- c(0,0,1,1,2,2)           + 301
  TUR.ctrl@obs.vars["NL_LPUE",ac(1)]            <- 0                        + 401
#  TUR.ctrl@obs.vars["IBTS_Q1",ac(1)]           <- 0                        + 501
  TUR.ctrl@cor.obs[]                            <- NA
#  TUR.ctrl@cor.obs["SNS",1:5]                   <- c(0,rep(1,4))
#  TUR.ctrl@cor.obs.Flag[2]                      <- af("AR")
  TUR.ctrl@biomassTreat[4]                      <- 2
  TUR.ctrl                                      <- update(TUR.ctrl)
  ### ------------------------------------------------------------------------------------------------------
  ###   4. Run assessment
  ### ------------------------------------------------------------------------------------------------------

  TUR.sam                 <- FLSAM(TUR,TUR.tun,TUR.ctrl)
  TUR.ctrl@residuals      <- FALSE; TUR.sam@control@residuals <- FALSE
  TUR.retro               <- retro(TUR,TUR.tun,TUR.ctrl,retro=5) # base.assess=TUR.sam
  TUR.sams[[ac(pg)]]      <- TUR.sam
  TUR.sams.retro[[ac(pg)]]<- TUR.retro
  TUR.ctrl@residuals      <- TRUE; TUR.sam@control@residuals <- TRUE
}

### ------------------------------------------------------------------------------------------------------
###   5. Diagnostics
### ------------------------------------------------------------------------------------------------------

for(pg in c(10:6)){
  sens <- paste0("pg",pg)
  TUR.sam <- TUR.sams[[ac(pg)]]
  TUR.retro <- TUR.sams.retro[[ac(pg)]]
  source(file.path(codePath,"03b_runDiagnostics.r"))
}

TUR.sams <- as(TUR.sams,"FLSAMs")
TUR.sams.retro <- as(TUR.sams.retro,"FLSAMs")

pdf(file.path(outPath,paste0(run,"_","pgcomb","assessmentOut.pdf")))
plot(TUR.sams)

par(mfrow=c(2,1))
print(plot(AIC(TUR.sams),ylab="AIC",xaxt="n",las=2,pch=19,xlab=""))
axis(1,at=1:5,labels=paste0("pg ",names(TUR.sams)),las=1)
print(grid())

print(plot(unlist(lapply(lapply(TUR.sams.retro,mohns.rho,ref.year=2017,span=5,type="ssb"),function(x){return(mean(x$rho[1:5]))})),xlab="",ylab="Mohns rho (5-year peel)",xaxt="n",las=2,pch=19,main="SSB"))
axis(1,at=1:5,labels=paste0("pg ",names(TUR.sams.retro)))
print(grid())

print(plot(unlist(lapply(lapply(TUR.sams.retro,mohns.rho,ref.year=2017,span=5,type="fbar"),function(x){return(mean(x$rho[1:5]))})),xlab="",ylab="Mohns rho (5-year peel)",xaxt="n",las=2,pch=19,main="Fbar"))
axis(1,at=1:5,labels=paste0("pg ",names(TUR.sams.retro)))
print(grid())

print(plot(unlist(lapply(lapply(TUR.sams.retro,mohns.rho,ref.year=2017,span=5,type="rec"),function(x){return(mean(x$rho[1:5]))})),xlab="",ylab="Mohns rho (5-year peel)",xaxt="n",las=2,pch=19,main="Recruitment"))
axis(1,at=1:5,labels=paste0("pg ",names(TUR.sams.retro)))
print(grid())

dev.off()

#- Compare mohns rho's of different runs
res <- rbind(unlist(lapply(lapply(TUR.sams.retro,mohns.rho,ref.year=2017,span=5,type="ssb"),function(x){return(mean(x$rho[1:5]))})),
             unlist(lapply(lapply(TUR.sams.retro,mohns.rho,ref.year=2017,span=5,type="fbar"),function(x){return(mean(x$rho[1:5]))})),
             unlist(lapply(lapply(TUR.sams.retro,mohns.rho,ref.year=2017,span=5,type="rec"),function(x){return(mean(x$rho[1:5]))})))
res           <- rbind(res,colMeans(abs(res)),colSums(abs(res)))
rownames(res) <- c("ssb","fbar","rec","mean","sum")
kable(res)

res <- rbind(unlist(lapply(lapply(TUR.sams.retro,mohns.rho,ref.year=2017,span=5,type="ssb"),function(x){return(mean(abs(x$rho[1:5])))})),
             unlist(lapply(lapply(TUR.sams.retro,mohns.rho,ref.year=2017,span=5,type="fbar"),function(x){return(mean(abs(x$rho[1:5])))})),
             unlist(lapply(lapply(TUR.sams.retro,mohns.rho,ref.year=2017,span=5,type="rec"),function(x){return(mean(abs(x$rho[1:5])))})))
res           <- rbind(res,colMeans(abs(res)),colSums(abs(res)))
rownames(res) <- c("ssb","fbar","rec","mean","sum")
kable(res)

#- Plot proportion of fish in plusgroup
storePG <- matrix(NA,nrow=length(6:10),ncol=length(1981:2017),dimnames=list(pg=6:10,year=1981:2017))
for(pg in 6:10)
  storePG[ac(pg),] <- sweep(TUR.sams[[ac(pg)]]@stock.n[ac(pg),],2:5,quantSums(TUR.sams[[ac(pg)]]@stock.n),"/")@.Data*100
matplot(y=t(storePG),x=1981:2017,type="b",pch=ac(c(6:9,0)),yaxs="i",ylab="Percentage of stock numbers in plusgroup",xlab="Years")
grid()