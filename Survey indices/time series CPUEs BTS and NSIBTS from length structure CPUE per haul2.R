ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}

specTranslate <- function(survdata){
  survdata$Species <- as.character(survdata$Species)
  try(survdata[survdata$Species %in% c("Dipturus batis","Dipturus flossada","Dipturus intermedia"), ]$Species <- "Dipturus batis complex", silent=T)

  #survdata[survdata$Species %in% c("Mustelus asterias","Mustelus mustelus","Mustelus"), ]$Species <- "Mustelus spp"

  #translate Sciliorhinus caniculus to Scyliorhinus canicula
  try(survdata[survdata$Species == "Sciliorhinus caniculus",]$Species <- "Scyliorhinus canicula", silent=T)

  #translate R brachyura to the two stocks
  try(survdata[survdata$Species == "Raja brachyura" & substr(survdata$SubArea,1,2) <= 36,]$Species <- "Raja brachyura 4c", silent=T)
  try(survdata[survdata$Species == "Raja brachyura" & substr(survdata$SubArea,1,2) >= 45,]$Species <- "Raja brachyura 4a", silent=T)
 
  try(survdata[survdata$Species %in% c("Psetta maxima","Scophthalmus maximus"),]$Species <- "Scophthalmus maximus", silent=T)
  
  return(survdata)
}


makeAnnualCPUE <- function(survdata,lwtable){
  
  #aggregate so that we aggregate complexes and sex
  survdata <- aggregate(CPUE_number_per_hour~Survey+ Ship + HaulNo + ShootLat +ShootLong + Year+Quarter+SubArea+Area+Species +LngtClass,data=survdata, sum)
 
  #make trawllist but make sure we do the right thing for the brachyura 4a and 4c later
  survdata$uniqueHaul <- paste(survdata$Year, survdata$Survey, survdata$Ship, survdata$HaulNo, survdata$ShootLat, survdata$ShootLong, sep="_")
  trawllistall <- survdata[!duplicated(survdata$uniqueHaul),!names( survdata) %in% c("Species","CPUE_number_per_hour","LngtClass")]
  
  
  CPUEsubarlengthAll <- CPUEsubarAll  <-   CPUEareaAll  <- CPUEyearAll <- NULL    

  #make selection
  for (spec in unique(survdata$Species)){
    cat(spec,"  ")
    datasel  <- survdata[survdata$Species==spec,]
    
    if(spec == "Raja brachyura 4a"){
      datasel <- datasel[substr(datasel$SubArea,1,2) >= 45,]
      trawllist <- trawllistall[substr(trawllistall$SubArea,1,2) >= 45,]
    }
    if(spec == "Raja brachyura 4c"){
      datasel <- datasel[substr(datasel$SubArea,1,2) <= 36,]
      trawllist <- trawllistall[substr(trawllistall$SubArea,1,2) <= 36,]
    } 
    if ( !(spec %in% c("Raja brachyura 4a","Raja brachyura 4c") ) ){
      trawllist <- trawllistall
    }
    
    if(exists("datasel")){
  
      # add zero catches for all lengths (adding shootlat and shootlon to unique haul is necessary because in 1985 IBTS Q1 Year,Survey,Ship,HaulNo is not unique )
      res <- expand.grid(uniqueHaul=unique(datasel$uniqueHaul),Species=unique(datasel$Species),LngtClass=unique(datasel$LngtClass), stringsAsFactors = F)
      cat("recs: ",nrow(res))
      res <- merge(datasel[,c("uniqueHaul","Species","LngtClass","CPUE_number_per_hour")],res, all=T)
      res <- merge(trawllist,res, all=T)
      try(res[is.na(res$CPUE_number_per_hour),]$CPUE_number_per_hour <- 0, silent=T)
  
      #add species (for missing hauls)
      res$Species <- spec
      try(res[is.na(res$LngtClass),]$LngtClass <- 0, silent =T)
  
      #aggregate to per rect length
      CPUEsubarlength <- aggregate(CPUE_number_per_hour~Survey+Year+Quarter+SubArea+Area+Species+LngtClass,data=res, mean)
      cat(", no recs CPUE per len rect",nrow(CPUEsubarlength))

      #merge with lwtable so that weights can also be calc'd
      CPUEsubarlength <- merge(CPUEsubarlength,lwtable)
      #a and b are units cm and gm, so divide lengths by 10 and res wts by 1000
      CPUEsubarlength$CPUE_kg_per_hour <- CPUEsubarlength$CPUE_number_per_hour * (CPUEsubarlength$a * (CPUEsubarlength$LngtClass/10)^CPUEsubarlength$b)/1000   
      #then estimate weights of exploitable biomass by taking only indvs over 30 cm (300 mm)
      CPUEsubarlength$exploitable <- 0
      if (any(CPUEsubarlength$LngtClass > 300)) CPUEsubarlength[CPUEsubarlength$LngtClass > 300,]$exploitable <- 1
      CPUEsubarlength$CPUE_kg_eb_per_hour <- CPUEsubarlength$exploitable * CPUEsubarlength$CPUE_kg_per_hour
      
      CPUEsubarlength$CPUE_kg_per_hour <- CPUEsubarlength$CPUE_number_per_hour * (CPUEsubarlength$a * (CPUEsubarlength$LngtClass/10)^CPUEsubarlength$b)/1000   
      
      
      #from per rect length to per rect
      CPUEsubar <- aggregate(cbind(CPUE_number_per_hour,CPUE_kg_per_hour,CPUE_kg_eb_per_hour)~Survey+Year+Quarter+SubArea+Area+Species,data=CPUEsubarlength, sum)
      cat(", no recs CPUE per rect",nrow(CPUEsubar))

      #from per rect to per area 
      CPUEarea <- aggregate(cbind(CPUE_number_per_hour,CPUE_kg_per_hour,CPUE_kg_eb_per_hour)~Survey+Year+Quarter+Area+Species,data=CPUEsubar, mean)
      cat(", no recs CPUE per area",nrow(CPUEarea))

      #from per area to per year, selecting only areas 1-7
      CPUEyear <- aggregate(cbind(CPUE_number_per_hour,CPUE_kg_per_hour,CPUE_kg_eb_per_hour)~Survey+Year+Quarter+Species,data=CPUEarea[CPUEarea$Area %in% 1:7,], mean)
      cat(", no recs CPUE per year",nrow(CPUEyear),"\n")
      if (is.null(CPUEyearAll)){CPUEyearAll <- CPUEyear}else{CPUEyearAll <- rbind(CPUEyearAll, CPUEyear)}
    }
  }
  return(CPUEyearAll)
}

maindir <- "w://imares/data/ices-wg/ibpturbot/2017/survey indices/"
#IBTS
survdataIBTS <- read.csv(paste0(maindir, "IBTS CPUE length haul files/CPUE per length per haul_2017-08-09 13_18_41.csv"))
#BTS
survdataBTS <- read.csv( paste0(maindir, "Datras CPUE BTS files/CPUE per length per Hour and Swept Area_2017-08-15 13_06_40.csv"))

#read LW table
lwtable <- read.csv(paste0(maindir,"lwtable.csv"))


survdataBTS$SubArea <- as.character(survdataBTS$StatRec)
survdataBTS$Area    <- 1                                                               #for BTS, act as if everything is in a single area that will be selected (because in IBTS roundfish areas 1-7 are selected)
# make good survey names
survdataBTS$Survey  <- paste(survdataBTS$Survey,survdataBTS$Country," ")
try(survdataBTS[survdataBTS$Country=="NED",]$Survey <- paste(survdataBTS[survdataBTS$Country=="NED",]$Survey,survdataBTS[survdataBTS$Country=="NED",]$Ship," "), silent=T)  #if lengthcass  is na that means it is zero
survdataBTS[is.na(survdataBTS$LngtClass),]$LngtClass <- 0
survdataBTS         <- survdataBTS[!(survdataBTS$HaulVal=="I"),]                               #remove invalid hauls
survdataBTS         <- survdataBTS[!( survdataBTS$Country=="ENG" & survdataBTS$Year < 1993) ,] #remove Uk data prior to 1993
survdataBTS         <- survdataBTS[survdataBTS$ICESArea %in% c("VIId","IVa","IVb","IVc"),] #only BTS data in IV and VIId (because ENG data started to have other areas in 2016)
survdataBTS         <- survdataBTS[!( survdataBTS$Country=="ENG" & survdataBTS$Year == 2016) ,] #remove Uk data in 2016 because CPUE calcs are off

#survdataCGFS$Area    <- 1 #for CGFS (like for BTS), act as if everything is in a single area that will be selected (because in IBTS roundfish areas 1-7 are selected)

survdataIBTS <- specTranslate(survdata=survdataIBTS)
survdataBTS  <- specTranslate(survdata=survdataBTS)
#survdataCGFS <- specTranslate(survdata=survdataCGFS)

toprint <- c("Scophthalmus maximus" )          

CPUEyearAllIBTS <- makeAnnualCPUE(survdata=survdataIBTS[survdataIBTS$Year %in% 1977:2017 & survdataIBTS$Species %in% toprint ,],lwtable)
CPUEyearAllIBTS$SQ <- paste(CPUEyearAllIBTS$Species, " ",CPUEyearAllIBTS$Survey," Q",CPUEyearAllIBTS$Quarter, sep="")
CPUEyearAllIBTS$CPUE <- CPUEyearAllIBTS$CPUE_number_per_hour

CPUEyearAllBTS  <- makeAnnualCPUE(survdata=survdataBTS[survdataBTS$Year %in% 1977:2017 & survdataBTS$Species %in% toprint ,],lwtable)
CPUEyearAllBTS$SQ <- paste(CPUEyearAllBTS$Species, " ",CPUEyearAllBTS$Survey," Q",CPUEyearAllBTS$Quarter, sep="")
CPUEyearAllBTS$CPUE <- CPUEyearAllBTS$CPUE_number_per_hour
CPUEyearAllBTS <- CPUEyearAllBTS[!CPUEyearAllBTS$Survey=="BTS BEL  ",] 

CPUEyearAllBTS <- CPUEyearAllBTS[CPUEyearAllBTS$Year > 1998,] #prior to 1998, there i for sure something wrong with the BTS ISIS indices in datras 



#CPUEyearAllCGFS <- makeAnnualCPUE(survdata=survdataCGFS[survdataCGFS$Year %in% 1977:2017 & survdataCGFS$Species %in% toprint ,],lwtable)
#CPUEyearAllCGFS$SQ <- paste(CPUEyearAllCGFS$Species, " ",CPUEyearAllCGFS$Survey," Q",CPUEyearAllCGFS$Quarter, sep="")
#CPUEyearAllCGFS$CPUE <- CPUEyearAllCGFS$CPUE_number_per_hour

CPUEyearAll <- rbind(CPUEyearAllIBTS,CPUEyearAllBTS) #, CPUEyearAllCGFS)

##############################################################
# Make tables
##############################################################

for (ii in c("CPUE_number_per_hour","CPUE_kg_per_hour","CPUE_kg_eb_per_hour")){
  w <- reshape(CPUEyearAll[,c("SQ","Year",ii)], 
             timevar = "SQ", idvar = c("Year"), direction = "wide")

  for (spec in toprint){
    print (spec);
    print(ii)
    repdat <-  w[,c(1,grep( spec, colnames(w), fixed=T))]
    names(repdat) <-  gsub(" ","",gsub(paste(ii,".",spec, sep=""),"",names(repdat)))
    write.csv(round(repdat,3), file=paste0(maindir,"IBTS output/",ii," ",spec,".csv"), row.names=F)
    print(round(repdat,3))
  }
}

##############################################################
# Make figures
##############################################################

lastyear <- 2017
#windows(9.5,12.5,record=T)
for(sp in unique(CPUEyearAllIBTS$Species)){
  #windows(12,21)
  par(mfrow=c(2,3),mar=c(3.0,3.5,0.3,0.5), mgp=c(2,0.8,0))
  cols= c('blue','red', 'black', 'purple')
  sqs <- unique(CPUEyearAllIBTS[CPUEyearAllIBTS$Species==sp,]$SQ)
  
  #first numbers
  
   if(sp=="Raja clavata"){ylim=c(0,6)}else{ylim=c(0,max(CPUEyearAllIBTS[CPUEyearAllIBTS$Species==sp,]$CPUE))}
   #ylim=c(0,1.1*max(CPUEyearAll[CPUEyearAll$Species==sp,]$CPUE))
  
   plotdat <- merge(CPUEyearAllIBTS[CPUEyearAllIBTS$Species==sp,],expand.grid(Year=1977:lastyear,SQ=sqs),by=c("Year", "SQ"), all=T) 
   plot(x=plotdat[plotdat$SQ==sqs[1],]$Year, y=plotdat[plotdat$SQ==sqs[1],]$CPUE, type="l", ylim=ylim,xlab="Year", ylab="CPUE (N per hour)", col=cols[1])
   lines(x=plotdat[plotdat$SQ==sqs[1],]$Year, y=ma(plotdat[plotdat$SQ==sqs[1],]$CPUE), lwd=3, col=cols[1])
   for (sq in sqs ){
     lines(x=plotdat[plotdat$SQ==sq,]$Year, y=plotdat[plotdat$SQ==sq,]$CPUE, type="l", ylim=ylim,xlab="Year", ylab="CPUE (N per hour)", col=cols[which(sq ==sqs)])
     lines(x=plotdat[plotdat$SQ==sq,]$Year, y=ma(plotdat[plotdat$SQ==sq,]$CPUE), lwd=3, col=cols[which(sq ==sqs)])
   }
   grid()
   legend("topright", sp, bty="n", inset=0.03, xjust=1)
   legend("topleft", unique(paste(plotdat[!is.na(plotdat$Survey),]$Survey," Q" ,plotdat[!is.na(plotdat$Survey),]$Quarter,sep="")), text.col=cols, bty="n", inset=0.03, xjust=0)

   #then IBTS wts
   if(sp=="Raja clavata"){ylim=c(0,6)}else{ylim=c(0,max(CPUEyearAllIBTS[CPUEyearAllIBTS$Species==sp,]$CPUE_kg_per_hour))}
   
   plotdat <- merge(CPUEyearAllIBTS[CPUEyearAllIBTS$Species==sp,],expand.grid(Year=1977:lastyear,SQ=sqs),by=c("Year", "SQ"), all=T) 
   plot(x=plotdat[plotdat$SQ==sqs[1],]$Year, y=plotdat[plotdat$SQ==sqs[1],]$CPUE_kg_per_hour, type="l", ylim=ylim,xlab="Year", ylab="CPUE (kg per hour)", col=cols[1])
   lines(x=plotdat[plotdat$SQ==sqs[1],]$Year, y=ma(plotdat[plotdat$SQ==sqs[1],]$CPUE_kg_per_hour), lwd=3, col=cols[1])
   for (sq in sqs ){
     lines(x=plotdat[plotdat$SQ==sq,]$Year, y=plotdat[plotdat$SQ==sq,]$CPUE_kg_per_hour, col=cols[which(sq ==sqs)])
     lines(x=plotdat[plotdat$SQ==sq,]$Year, y=ma(plotdat[plotdat$SQ==sq,]$CPUE_kg_per_hour), lwd=3, col=cols[which(sq ==sqs)])
   }
   grid()
   legend("topright", sp, bty="n", inset=0.03)
   legend("topleft", unique(paste(plotdat[!is.na(plotdat$Survey),]$Survey," Q" ,plotdat[!is.na(plotdat$Survey),]$Quarter,sep="")), text.col=cols, bty="n", inset=0.03)
   
   #then IBTS exploitable biomass wts
   if(sp=="Raja clavata"){ylim=c(0,6)}else{ylim=c(0,max(CPUEyearAllIBTS[CPUEyearAllIBTS$Species==sp,]$CPUE_kg_eb_per_hour))}
   
   plotdat <- merge(CPUEyearAllIBTS[CPUEyearAllIBTS$Species==sp,],expand.grid(Year=1977:lastyear,SQ=sqs),by=c("Year", "SQ"), all=T) 
   plot(x=plotdat[plotdat$SQ==sqs[1],]$Year, y=plotdat[plotdat$SQ==sqs[1],]$CPUE_kg_eb_per_hour, type="l", ylim=ylim,xlab="Year", ylab="CPUE > 30 cm (kg per hour)", col=cols[1])
   lines(x=plotdat[plotdat$SQ==sqs[1],]$Year, y=ma(plotdat[plotdat$SQ==sqs[1],]$CPUE_kg_eb_per_hour), lwd=3, col=cols[1])
   for (sq in sqs ){
     lines(x=plotdat[plotdat$SQ==sq,]$Year, y=plotdat[plotdat$SQ==sq,]$CPUE_kg_eb_per_hour, col=cols[which(sq ==sqs)])
     lines(x=plotdat[plotdat$SQ==sq,]$Year, y=ma(plotdat[plotdat$SQ==sq,]$CPUE_kg_eb_per_hour), lwd=3, col=cols[which(sq ==sqs)])
   }
   grid()
   legend("topright", sp, bty="n", inset=0.03)
   legend("topleft", unique(paste(plotdat[!is.na(plotdat$Survey),]$Survey," Q" ,plotdat[!is.na(plotdat$Survey),]$Quarter,sep="")), text.col=cols, bty="n", inset=0.03)
   
   
   #then BTS numbers
  ylim=c(0,1.1*max(CPUEyearAllBTS[CPUEyearAllBTS$Species==sp,]$CPUE))
  sqs <- unique(CPUEyearAllBTS[CPUEyearAllBTS$Species==sp,]$SQ)
  
  plotdat <- merge(CPUEyearAllBTS[CPUEyearAllBTS$Species==sp,],expand.grid(Year=1977:lastyear,SQ=sqs),by=c("Year", "SQ"), all=T) 
  if (nrow(plotdat)>0){
    plot(x=plotdat[plotdat$SQ==sqs[1],]$Year, y=plotdat[plotdat$SQ==sqs[1],]$CPUE, type="l", ylim=ylim,xlab="Year", ylab="CPUE (N per hour)", col=cols[1])
    lines(x=plotdat[plotdat$SQ==sqs[1],]$Year, y=ma(plotdat[plotdat$SQ==sqs[1],]$CPUE), lwd=3, col=cols[1])
    for (sq in sqs ){
      lines(x=plotdat[plotdat$SQ==sq,]$Year, y=plotdat[plotdat$SQ==sq,]$CPUE, type="l", col=cols[which(sq ==sqs)])
      lines(x=plotdat[plotdat$SQ==sq,]$Year, y=ma(plotdat[plotdat$SQ==sq,]$CPUE), lwd=3, col=cols[which(sq ==sqs)])
    }
    grid()
    legend("topright", sp, bty="n", inset=0.03)
    legend("topleft", unique(paste(plotdat[!is.na(plotdat$Survey),]$Survey," Q" ,plotdat[!is.na(plotdat$Survey),]$Quarter,sep="")), text.col=cols, bty="n", inset=0.03)
  }
  
  #then bts wts
  ylim=c(0,1.1*max(CPUEyearAllBTS[CPUEyearAllBTS$Species==sp,]$CPUE_kg_per_hour))
  
  plotdat <- merge(CPUEyearAllBTS[CPUEyearAllBTS$Species==sp,],expand.grid(Year=1977:lastyear,SQ=sqs),by=c("Year", "SQ"), all=T) 
  if (nrow(plotdat)>0){
    plot(x=plotdat[plotdat$SQ==sqs[1],]$Year, y=plotdat[plotdat$SQ==sqs[1],]$CPUE_kg_per_hour, type="l", ylim=ylim,xlab="Year", ylab="CPUE (kg per hour)", col=cols[1])
    lines(x=plotdat[plotdat$SQ==sqs[1],]$Year, y=ma(plotdat[plotdat$SQ==sqs[1],]$CPUE_kg_per_hour), lwd=3, col=cols[1])
    for (sq in sqs ){
      lines(x=plotdat[plotdat$SQ==sq,]$Year, y=plotdat[plotdat$SQ==sq,]$CPUE_kg_per_hour, type="l", col=cols[which(sq ==sqs)])
      lines(x=plotdat[plotdat$SQ==sq,]$Year, y=ma(plotdat[plotdat$SQ==sq,]$CPUE_kg_per_hour), lwd=3, col=cols[which(sq ==sqs)])
    }
    grid()
    legend("topright", sp, bty="n", inset=0.03)
    legend("topleft", unique(paste(plotdat[!is.na(plotdat$Survey),]$Survey," Q" ,plotdat[!is.na(plotdat$Survey),]$Quarter,sep="")), text.col=cols, bty="n", inset=0.03)
  }
 
   #then bts exp biomass wts
  ylim=c(0,1.1*max(CPUEyearAllBTS[CPUEyearAllBTS$Species==sp,]$CPUE_kg_eb_per_hour))
  
  plotdat <- merge(CPUEyearAllBTS[CPUEyearAllBTS$Species==sp,],expand.grid(Year=1977:lastyear,SQ=sqs),by=c("Year", "SQ"), all=T) 
  if (nrow(plotdat)>0){
    plot(x=plotdat[plotdat$SQ==sqs[1],]$Year, y=plotdat[plotdat$SQ==sqs[1],]$CPUE_kg_eb_per_hour, type="l", ylim=ylim,xlab="Year", ylab="CPUE > 50 cm (kg per hour)", col=cols[1])
    lines(x=plotdat[plotdat$SQ==sqs[1],]$Year, y=ma(plotdat[plotdat$SQ==sqs[1],]$CPUE_kg_eb_per_hour), lwd=3, col=cols[1])
    for (sq in sqs ){
      lines(x=plotdat[plotdat$SQ==sq,]$Year, y=plotdat[plotdat$SQ==sq,]$CPUE_kg_eb_per_hour, type="l", col=cols[which(sq ==sqs)])
      lines(x=plotdat[plotdat$SQ==sq,]$Year, y=ma(plotdat[plotdat$SQ==sq,]$CPUE_kg_eb_per_hour), lwd=3, col=cols[which(sq ==sqs)])
    }
    grid()
    legend("topright", sp, bty="n", inset=0.03)
    legend("topleft", unique(paste(plotdat[!is.na(plotdat$Survey),]$Survey," Q" ,plotdat[!is.na(plotdat$Survey),]$Quarter,sep="")), text.col=cols, bty="n", inset=0.03)
  }
  
  
  #then cgfs nums
#  ylim=c(0,1.1*max(CPUEyearAllCGFS[CPUEyearAllCGFS$Species==sp,]$CPUE))
  
#  sqs <- unique(CPUEyearAllCGFS[CPUEyearAllCGFS$Species==sp,]$SQ)
#  cols= c('blue','red', 'black', 'purple')
#  plotdat <- merge(CPUEyearAllCGFS[CPUEyearAllCGFS$Species==sp,],expand.grid(Year=1977:lastyear,SQ=sqs),by=c("Year", "SQ"), all=T) 
#  if (nrow(plotdat)>0){
#    plot(x=plotdat[plotdat$SQ==sqs[1],]$Year, y=plotdat[plotdat$SQ==sqs[1],]$CPUE, type="l", ylim=ylim,xlab="Year", ylab="CPUE (N per hour)", col=cols[1])
#    lines(x=plotdat[plotdat$SQ==sqs[1],]$Year, y=ma(plotdat[plotdat$SQ==sqs[1],]$CPUE), lwd=3, col=cols[1])
#    for (sq in sqs ){
#      lines(x=plotdat[plotdat$SQ==sq,]$Year, y=plotdat[plotdat$SQ==sq,]$CPUE, type="l", col=cols[which(sq ==sqs)])
#      lines(x=plotdat[plotdat$SQ==sq,]$Year, y=ma(plotdat[plotdat$SQ==sq,]$CPUE), lwd=3, col=cols[which(sq ==sqs)])
#    }
#    grid()
#    legend("topright", sp, bty="n", inset=0.03)
#    legend("topleft", unique(paste(plotdat[!is.na(plotdat$Survey),]$Survey," Q" ,plotdat[!is.na(plotdat$Survey),]$Quarter,sep="")), text.col=cols, bty="n", inset=0.03)

    #then cgfs wts
#    ylim=c(0,1.1*max(CPUEyearAllCGFS[CPUEyearAllCGFS$Species==sp,]$CPUE_kg_per_hour))
    
#    sqs <- unique(CPUEyearAllCGFS[CPUEyearAllCGFS$Species==sp,]$SQ)
#    cols= c('blue','red', 'black', 'purple')
#    plotdat <- merge(CPUEyearAllCGFS[CPUEyearAllCGFS$Species==sp,],expand.grid(Year=1977:lastyear,SQ=sqs),by=c("Year", "SQ"), all=T) 
#    if (nrow(plotdat)>0){
#      plot(x=plotdat[plotdat$SQ==sqs[1],]$Year, y=plotdat[plotdat$SQ==sqs[1],]$CPUE_kg_per_hour, type="l", ylim=ylim,xlab="Year", ylab="CPUE > 50 cm (kg per hour)", col=cols[1])
#      lines(x=plotdat[plotdat$SQ==sqs[1],]$Year, y=ma(plotdat[plotdat$SQ==sqs[1],]$CPUE_kg_per_hour), lwd=3, col=cols[1])
#      for (sq in sqs ){
#        lines(x=plotdat[plotdat$SQ==sq,]$Year, y=plotdat[plotdat$SQ==sq,]$CPUE_kg_per_hour, type="l", col=cols[which(sq ==sqs)])
#        lines(x=plotdat[plotdat$SQ==sq,]$Year, y=ma(plotdat[plotdat$SQ==sq,]$CPUE_kg_per_hour), lwd=3, col=cols[which(sq ==sqs)])
#      }
#      grid()
#      legend("topright", sp, bty="n", inset=0.06)
#      legend("topleft", unique(paste(plotdat[!is.na(plotdat$Survey),]$Survey," Q" ,plotdat[!is.na(plotdat$Survey),]$Quarter,sep="")), text.col=cols, bty="n", inset=0.06)
#    }
    #then cgfs exploitable wts
#    ylim=c(0,1.1*max(CPUEyearAllCGFS[CPUEyearAllCGFS$Species==sp,]$CPUE_kg_eb_per_hour))
    
#    sqs <- unique(CPUEyearAllCGFS[CPUEyearAllCGFS$Species==sp,]$SQ)
#    cols= c('blue','red', 'black', 'purple')
#    plotdat <- merge(CPUEyearAllCGFS[CPUEyearAllCGFS$Species==sp,],expand.grid(Year=1977:lastyear,SQ=sqs),by=c("Year", "SQ"), all=T) 
#    if (nrow(plotdat)>0){
#      plot(x=plotdat[plotdat$SQ==sqs[1],]$Year, y=plotdat[plotdat$SQ==sqs[1],]$CPUE_kg_eb_per_hour, type="l", ylim=ylim,xlab="Year", ylab="CPUE (kg per hour)", col=cols[1])
#      lines(x=plotdat[plotdat$SQ==sqs[1],]$Year, y=ma(plotdat[plotdat$SQ==sqs[1],]$CPUE_kg_eb_per_hour), lwd=3, col=cols[1])
#      for (sq in sqs ){
#        lines(x=plotdat[plotdat$SQ==sq,]$Year, y=plotdat[plotdat$SQ==sq,]$CPUE_kg_eb_per_hour, type="l", col=cols[which(sq ==sqs)])
#        lines(x=plotdat[plotdat$SQ==sq,]$Year, y=ma(plotdat[plotdat$SQ==sq,]$CPUE_kg_eb_per_hour), lwd=3, col=cols[which(sq ==sqs)])
#      }
#      grid()
#      legend("topright", sp, bty="n", inset=0.03)
#      legend("topleft", unique(paste(plotdat[!is.na(plotdat$Survey),]$Survey," Q" ,plotdat[!is.na(plotdat$Survey),]$Quarter,sep="")), text.col=cols, bty="n", inset=0.03)
#    }
#  }
}
