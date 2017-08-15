#-------------------------------------------------------------------------------
#
# Disclaimer:
# Nobody is allowed to use this data without prior permission from CVO and a
# signed decleration on the usage of the VISSTAT database.
# Ask Peter / Sieto for guidance.
#
# Although I (Niels) am rather confident of the quality of this data, everyone
# is responsible for a thorough screening of the data. Please report errors to me.
#-------------------------------------------------------------------------------

#- Clear workspace
rm(list=ls())

library(vmstools) #- code.google.com/p/vmstools
library(maps); library(mapdata); library(rgdal); library(Matrix); library(data.table)

#- Define relevant paths
if(!substr(R.Version()$os,1,3)== "lin") sysPa<-paste("w",":/IMARES/data/",sep="")
if( substr(R.Version()$os,1,3)== "lin") sysPa<-"~/wur/N/Projecten/"
ProjectPa<-paste(sysPa,"vms_/",sep="")
codePath  <- paste(ProjectPa,"BASEFILES/",sep="")
dataPath  <- paste(ProjectPa,"BASEFILES/",sep="")
outPath   <- paste(ProjectPa,"BASEFILES/",sep="")
outPath   <- "w:/IMARES/Data/ICES-WG/IBPTURBOT/2017/LPUE/CatchEffort/Output/"
codePath  <- outPath 

#-------------------------------------------------------------------------------
#- 1) Load the data
#-------------------------------------------------------------------------------
load(file=file.path(outPath,paste("eflaloTTY9516.RData",sep="")));
eflaloTTY$year<-an(substr(eflaloTTY$LE_CDATIM,1,4))
eflaloTTY<-eflaloTTY[eflaloTTY$KwDays>0,]
eflaloTTY<-as.data.frame(cbind(eflaloTTY,ICESrectangle2LonLat(eflaloTTY$LE_RECT,midpoint=T)))
eflaloTTY$VS<-"large"
eflaloTTY[eflaloTTY$VE_KW<222,]$VS<-"Euro"
eflaloTTY<-eflaloTTY[!eflaloTTY$LE_GEAR_INNOV=="TBS_",]
eflaloTTY[grepl('PU', eflaloTTY$LE_GEAR_INNOV),]$LE_GEAR_INNOV<-"PULS"
eflaloTTY[grepl('SUM', eflaloTTY$LE_GEAR_INNOV),]$LE_GEAR_INNOV<-"SUM"
eflaloTTY[grepl('slof', eflaloTTY$LE_GEAR_INNOV),]$LE_GEAR_INNOV<-"SLOF"
unique(eflaloTTY$LE_GEAR_INNOV)

aggdata<-with(eflaloTTY, aggregate(list(LE_KG_TUR,LE_KG_SOL,LE_KG_PLE,KwDays,TRIP),
                                   by=list(LE_GEAR_INNOV,VS,year,SI_LATI, SI_LONG,LE_RECT), FUN=sum,na.rm=TRUE))
colnames(aggdata)<-c("LE_GEAR_INNOV","VS","year","SI_LATI","SI_LONG","LE_RECT","LE_KG_TUR","LE_KG_SOL","LE_KG_PLE","KwDays","TRIP")

unique(aggdata$LE_GEAR_INNOV)
head(aggdata)
save(aggdata,file=file.path(outPath,paste("aggdata9516.RData",sep="")))

aggdata$cols<-paste(aggdata$VS,aggdata$LE_GEAR_INNOV,sep="_")
tabel<-as.data.frame(with(aggdata, tapply(TRIP,list(year,cols),sum)))
tabel[is.na(tabel)]<-0
tabel$ALL<-rowSums(tabel)
tabel<-rbind(tabel,round(colMeans(tabel)))
row.names(tabel)[17]<-"Mean"
tabel.trips<-tabel
tabel<-as.data.frame(with(aggdata, tapply(KwDays,list(year,cols),sum)))
tabel[is.na(tabel)]<-0
tabel$ALL<-rowSums(tabel)
tabel<-rbind(tabel,round(colMeans(tabel)))
row.names(tabel)[17]<-"Mean"
tabel.KwDays<-tabel
tabel.KwDays<-round(tabel.KwDays,0)
tabel<-as.data.frame(with(aggdata, tapply(LE_KG_TUR,list(year,cols),sum)))
tabel[is.na(tabel)]<-0
tabel$ALL<-rowSums(tabel)
tabel<-rbind(tabel,round(colMeans(tabel)))
row.names(tabel)[17]<-"Mean"
tabel.LE_KG_TUR<-tabel
tabel.LpUE<-round(tabel.LE_KG_TUR/tabel.KwDays,3)
save(tabel.trips,tabel.KwDays,tabel.LE_KG_TUR,tabel.LpUE,file=file.path(outPath,paste("Tabels9516.RData",sep="")))
