#install.packages('DATRAS',repos='http://www.rforge.net/',type='source')
#library(devtools)
#install_github("casperwberg/surveyIndex/surveyIndex")


library(DATRAS)
library(mgcv);
library(parallel);
library(maps); library(mapdata);
library(surveyIndex);
## Species specific parameters:
cmSize=1;
spectrumMax=70;
agesQ1=1:10
agesQ3=1:10
years=1991:2016 
outFolder=".";
genus="Psetta"
bfamily="maxima";
path <- "w://imares/data/ICES-WG/IBPturbot_2017/survey indices"
setwd(path)

dAll <- readExchangeDir(paste0(path,"/exchange BTS/"),strict=FALSE)
#dAll<-addSpatialData(dAll,"../../shapefiles/ICES_areas.shp")
dAll<-subset(dAll,Species==paste(genus,bfamily),Year %in% years,HaulVal=="V",StdSpecRecCode==1)
  
dAll=addSpectrum(dAll,cm.breaks=seq(0,spectrumMax,by=cmSize))

## impute missing depths
summary(dAll$Depth)

dmodel=gam(log(Depth) ~ s(lon,lat,k=200),data=dAll[[2]])
sel=subset(dAll,is.na(Depth))
sel$Depth=0; ## Guard against NA-error
dAll$Depth[is.na(dAll$Depth)]=exp(predict(dmodel,newdata=sel[[2]]))
dmodel=NULL
sel=NULL
gc()


dQ3=dAll

dQ3=addWeightByHaul(dQ3)

dQ3.BTS = subset(dQ3,Survey=="BTS")

mybubblePlot<-function (d, response = "HaulWgt", scale = NULL, col.zero = "red", 
                        pch.zero = "+", ...) 
{
  d[[2]]$resp.var <- d[[2]][[response]]
  if (is.null(scale)) 
    scale = mean(d[[2]]$resp.var, na.rm = TRUE)/max(d[[2]]$resp.var, 
                                                    na.rm = TRUE)
  plot(d$lon, d$lat, type = "n", xlab = "Longitude", ylab = "Latitude",...)
  map("worldHires", fill = TRUE, plot = TRUE, add = TRUE, col = grey(0.5))
  points(d$lon, d$lat, pch = 16, cex = scale * sqrt(d[[2]]$resp.var), 
         ...)
  zero = subset(d, resp.var == 0)
  points(zero$lon, zero$lat, pch = pch.zero, col = col.zero)
}

par(mfrow=c(1,3))
mybubblePlot(dQ3.BTS,scale=1/5,ylim=c(51,62),main="BTS Q3")
rect(-0.7,58.24,5,62,border=4,lwd=2)

## Area subsetting
dQ3=subset(dQ3,!(lon>-0.7 & lon<5 & lat>58.24 & lat<62))
dQ3.BTS=subset(dQ3.BTS,!(lon>-0.7 & lon<5 & lat>58.24 & lat<62))

## tables
xtabs(~Year+Gear,data=dQ3[[2]])
xtabs(~Year+Country,data=dQ3[[2]])


library(xtable)
sink("tables.txt")
print(xtable(xtabs(~Year+Gear,data=dQ3.BTS[[2]]),digits=0,caption="Number of hauls by gear and year BTS Q3"))

## Gear subsetting
dQ3=subset(dQ3,!Gear %in% c("BT4S","ABD"))
dQ3.BTS=subset(dQ3.BTS,!Gear %in% c("BT4S","ABD"))

print(xtable(xtabs(NoAtALK~Year+Survey,data=dQ3[[1]]),digits=0,caption="Number of age samples by year and survey (excluding BT4S and ABD gears) Q3"))

print(xtable(xtabs(NoAtALK~Year+Age,data=dQ3.BTS[[1]]),digits=0,caption="Number of age samples by year and age in BTS (excluding BT4S gear) Q3"))

print(xtable(xtabs(NoAtALK~Year+Age,data=subset(dQ3.BTS[[1]],Age<13)),digits=0,caption="Number of age samples by year and age in BTS (excluding BT4S gear) Q3"))

print(xtable(xtabs(NoAtALK~Year+Age,data=subset(dQ3.IBTS[[1]],Age<13)),digits=0,caption="Number of age samples by year and age in IBTS (excluding ABD gear) Q3"))

tmp=merge(dQ3.IBTS[[1]],dQ3.IBTS[[2]][,c("haul.id","ICESAREA")],by="haul.id",all.x=TRUE)

print(xtable(xtabs(NoAtALK~Year+ICESAREA,data=tmp),digits=0,caption="Number of age samples by year and ICES area in IBTS (excluding ABD gear) Q3"))

sink()

tmp=NULL
gc()

## Year subsetting
dQ3=subset(dQ3,Year %in% 1996:2016)
dQ3.BTS=subset(dQ3.BTS,Year %in% 1996:2016)
dQ3.IBTS=subset(dQ3.IBTS,Year %in% 1996:2016)

removeAgeNAs<-function(x) {
  x[[1]]=subset(x[[1]],!is.na(x[[1]]$Age))
  x[[1]]=subset(x[[1]],!is.na(x[[1]]$NoAtALK))
  x
}

dQ3=removeAgeNAs(dQ3)
dQ3.BTS=removeAgeNAs(dQ3.BTS)
dQ3.IBTS=removeAgeNAs(dQ3.IBTS)

##################

###############
## Declare settings for ALK model
mf = "" 
ack=TRUE;
useBICs=TRUE;
varCofs=FALSE;
maxKs=50;
mc.cores=1

add.ALK<-function(d){
  
  ages=agesQ3
  
  if(d$Quarter[1]=="1"){
    d[[1]]=subset(d[[1]],Age>0)
    d=fixAgeGroup(d,1)
    ages=agesQ1
  }
  
  d=addSpectrum(d,cm.breaks=seq(0,spectrumMax,by=cmSize))
  
  d.ysplit = split(d,d$Year)
  
  d.ALK= mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=ack,useBIC=useBICs,varCof=varCofs,maxK=maxKs,mc.cores=mc.cores)
  
  d.Nage=mclapply(d.ALK,predict,mc.cores=mc.cores)
  for(i in 1:length(d.ALK)) d.ysplit[[i]]$Nage=d.Nage[[i]];
  dd <- do.call("c",d.ysplit)
  dd    
}

dQ3=add.ALK(dQ3)
dQ3.BTS=add.ALK(dQ3.BTS)

xtabs(NoAtALK ~ Year+Age,data=dQ3.IBTS[[1]])

grid.BTS=getGrid(dQ3.BTS,nLon=40) 
grid.comb=getGrid(dQ3,nLon=40)

BTSmodels=list()
combmodels=list()
IBTSmodels=list()


## Make ctime : a numeric time variable 
dQ3$ctime = as.numeric(as.character(dQ3$Year))
dQ3.BTS$ctime = as.numeric(as.character(dQ3.BTS$Year))


#######################
## Model formulae
#######################

## Stationary model
modelsStatZ=rep("Year+Gear+s(lon,lat,bs=c('tp'),k=kvecZ[a])+s(Depth,bs='ts',k=6)+offset(log(HaulDur))",length(agesQ3))

modelsStatP=rep("Year+Gear+s(lon,lat,bs=c('tp'),k=kvecP[a])+s(Depth,bs='ts',k=6)+offset(log(HaulDur))",length(agesQ3))

modelsNonStat=rep("Year+Gear+te(ctime,lon,lat,d=c(1,2),bs=c('cs','tp'),k=c(5,25))+s(Depth,bs='ts',k=6)+offset(log(HaulDur))",length(agesQ3))



mc.cores=1
BTSmodels$SI.ac=getSurveyIdx(dQ3.BTS,ages=0:10,myids=grid.BTS[[3]],cutOff=0.5,fam="LogNormal",mc.cores=mc.cores,modelZ=modelsStatZ,modelP=modelsStatP)

exportSI<-function(x,ages,years,toy,file,nam="",exclude=c()){
  cat(nam,"\n",file=file)
  cat(range(as.numeric(as.character(years))),"\n",file=file,append=TRUE)
  cat("1 1 ",rep(toy,2),"\n",file=file,append=TRUE)
  cat(min(ages),max(ages),"\n",file=file,append=TRUE)
  write.table(round(cbind(1,x[,]),4),file=file,row.names=FALSE,col.names=FALSE,append=TRUE)
}

exportSI(BTSmodels$SI.ac$idx,agesQ3,1996:2016,toy=mean(dQ3.BTS[[2]]$timeOfYear,na.rm=TRUE),file="Q3BTSx2.dat",nam=paste("NS Plaice ; Last age is plus group, calculated",Sys.time())); 


summary(dQ3.BTS[[2]]$Country)


## Remove Belgian data - very little impact
##dQ3.BTS=subset(dQ3.BTS,Country!="BEL")

##BTSmodels$SI.ac.noBEL=getSurveyIdx(dQ3.BTS,ages=agesQ3,myids=grid.BTS[[3]],cutOff=0.5,fam="LogNormal",mc.cores=mc.cores,modelZ=modelsStatZ,modelP=modelsStatP)

##BTSmodels$SI.ac.noBEL<-NULL
##gc()



####################
## Combined
####################

combmodels$SI.ac=getSurveyIdx(dQ3,ages=agesQ3,myids=grid.comb[[3]],cutOff=0.5,fam="LogNormal",mc.cores=2,modelZ=modelsStatZ,modelP=modelsStatP)

exportSI(combmodels$SI.ac$idx,agesQ3,1996:2016,toy=mean(dQ3[[2]]$timeOfYear,na.rm=TRUE),file="Q3comb.dat",nam=paste("NS Plaice ; Last age is plus group, calculated",Sys.time())); 


#############################
## IBTS + BTS age samples
#############################
## only keep age samples from BTS in IBTS Q3 series
dQ3.IBTS=dQ3
dQ3.IBTS[[2]]<-subset(dQ3.IBTS[[2]],Survey=="NS-IBTS")
dQ3.IBTS[[3]]<-subset(dQ3.IBTS[[3]],Survey=="NS-IBTS")

grid.IBTS=getGrid(dQ3.IBTS,nLon=40)

IBTSmodels$SI.ac=getSurveyIdx(dQ3.IBTS,ages=0:10,myids=grid.IBTS[[3]],cutOff=0.15,fam="LogNormal",mc.cores=mc.cores,modelZ=modelsStatZ,modelP=modelsStatP)

exportSI(IBTSmodels$SI.ac$idx,agesQ3,1996:2016,toy=mean(dQ3.IBTS[[2]]$timeOfYear,na.rm=TRUE),file="Q3IBTS.dat",nam=paste("NS Plaice ; Last age is plus group, calculated",Sys.time())); 
