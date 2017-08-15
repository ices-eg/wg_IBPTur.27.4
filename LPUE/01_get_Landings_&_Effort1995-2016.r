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


for(year in ac(1995:2016)){
  print(year)
  #-------------------------------------------------------------------------------
  #- 1c) load tacsat and eflalo data from file
  #-------------------------------------------------------------------------------
  if(year >= 2001)
    load(file.path(dataPath,paste("cleanEflaloInnov",year,".RData",sep="")));
  if(year < 2001){
    load(file.path(dataPath,paste("cleanEflalo",year,".RData",sep="")));
    eflalo$LE_CDATIM  <- as.POSIXct(eflalo$LE_CDAT,format="%d/%m/%Y")
    eflalo$ELID       <- paste(eflalo$VE_COU,eflalo$VE_REF,eflalo$LE_GEAR,eflalo$LE_RECT)
    eflalo$LE_MONTH   <- format(eflalo$LE_CDATIM,"%m")
    eflalo$INTVF      <- eflalo$INTV
    eflalo$LE_GEAR_INNOV <- "TBB"
  }
  #head(eflalo[c("FT_DDATIM","FT_LDATIM","ELID","ID","LE_MONTH","INTV","INTVF","LE_GEAR_INNOV","LE_KG_TUR","LE_KG_SOL","LE_KG_PLE")])
  
  eflaloT<-cbind(eflalo[c("FT_REF",        "VE_REF",        "LE_CDATIM",     "VE_ID","LE_DIV",
                      "VE_COU",        "VE_LEN",        "VE_KW",         "VE_TON",
                      "FT_DCOU",       "FT_DHAR",       "FT_DDAT",       "FT_DTIME",
                      "FT_LCOU",       "FT_LHAR",       "FT_LDAT",       "FT_LTIME",
                      "LE_ID",         "LE_CDAT",       "LE_STIME",      "LE_ETIME",
                      "LE_SLAT",       "LE_SLON",       "LE_ELAT",       "LE_ELON",
                      "LE_GEAR",       "LE_WIDTH",      "LE_MSZ",        "LE_RECT")],
               eflalo[c("FT_DDATIM","FT_LDATIM","ELID","ID","LE_MONTH","INTV","INTVF","LE_GEAR_INNOV","LE_KG_TUR","LE_KG_SOL","LE_KG_PLE")])
  eflaloT$LE_MSZ<-an(eflaloT$LE_MSZ)
  eflaloT<-eflaloT[eflaloT$LE_GEAR=="TBB"&an(eflaloT$LE_MSZ)>70&an(eflaloT$LE_MSZ)<100&eflaloT$LE_DIV=="IV"&!is.na(eflaloT$ID),]
  eflaloT$LE_KG_TUR[is.na(eflaloT$LE_KG_TUR)]<-0
  eflaloT$LE_KG_SOL[is.na(eflaloT$LE_KG_SOL)]<-0
  eflaloT$LE_KG_PLE[is.na(eflaloT$LE_KG_PLE)]<-0
  eflaloT$INTV       <- an(difftime(eflaloT$FT_LDATIM,eflaloT$FT_DDATIM,units="days"))
  eflaloT$dummy      <- 1
  res               <- aggregate(eflaloT$dummy,by=list(eflaloT$FT_REF),FUN=sum,na.rm=T)
  colnames(res)     <- c("FT_REF","nrRecords")
  eflaloT            <- merge(eflaloT,res,by="FT_REF")
  eflaloT$INTV       <- eflaloT$INTV / eflaloT$nrRecords
  eflaloT$KwDays     <- eflaloT$INTV *eflaloT$VE_KW
  eflaloT$TRIP       <- 1 / eflaloT$nrRecords
  if(length(which(is.na(eflaloT$LE_GEAR_INNOV)))>0)
    eflaloT[is.na(eflaloT$LE_GEAR_INNOV),]$LE_GEAR_INNOV<-"TBB"
  eflaloTT<-eflaloT
  if (year==1995) eflaloTTY<-eflaloTT
                ColNames<-colnames(eflaloTTY)
  if (!(year==1995)) eflaloTTY<-rbind(eflaloTTY,eflaloTT[ColNames])
}
save(eflaloTTY,file=file.path(outPath,paste("eflaloTTY9516.RData",sep="")));
