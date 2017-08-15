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
library(RColorBrewer);library(ggplot2); library(lattice)
library(mgcv)
#- Define relevant paths
#- Define relevant paths
if(!substr(R.Version()$os,1,3)== "lin") sysPa<-paste("w",":/IMARES/data/",sep="")
if( substr(R.Version()$os,1,3)== "lin") sysPa<-"~/wur/N/Projecten/"
ProjectPa<-paste(sysPa,"vms_/",sep="")
codePath  <- paste(ProjectPa,"BASEFILES/",sep="")
dataPath  <- paste(ProjectPa,"BASEFILES/",sep="")
outPath   <- paste(ProjectPa,"BASEFILES/",sep="")
outPath   <- "w:/IMARES/Data/ICES-WG/IBPTURBOT/2017/LPUE/CatchEffort/Output/"
codePath  <- outPath 
load(file=file.path(outPath,paste("aggdata9516.RData",sep="")));
agglarge<-aggdata[aggdata$VS=="large",]
agglarge$LpUE<-agglarge$LE_KG_TUR/agglarge$KwDays
agglarge$value<-log(agglarge$LpUE+exp(-5))
agglarge<-agglarge[agglarge$SI_LATI<58,]
grs<-unique(agglarge$LE_GEAR_INNOV)


#- Make a figure of the LPUE by gear type
plotdat <- aggregate(agglarge[,"LpUE"],by=as.list(agglarge[,c("year","LE_GEAR_INNOV")]),FUN=median,na.rm=T)
xyplot(x ~ year,data=plotdat,group=LE_GEAR_INNOV,lty=1,type="l",lwd=2,auto.key=list(type="l",space = "right",lty=2,lwd=2),ylab="LpUE")


head(agglarge)

####
# remove 31F6 and 36F8 rectangle (GErmany, plaice box)
#########
agglarge <- agglarge[!agglarge$LE_RECT %in% c("31F6","36F8","36F7","36F8" ,"37F8", "35F7"),]

########################
# which years are sufficiently sampled
########################
barplot(sort(apply((table(agglarge$year,agglarge$LE_RECT))==0,2,sum)))
apply((table(agglarge$year,agglarge$LE_RECT))==0,2,sum)
lt4misdata <- names(which((apply((table(agglarge$year,agglarge$LE_RECT))==0,2,sum)) < 4))


gam1a       <- gam(value ~ te(SI_LONG,SI_LATI,k=5)+as.factor(year) + LE_GEAR_INNOV,weights=TRIP, data=agglarge); 
gam.check(gam1a)
summary(gam1a)
plot(gam1a, all.terms=T)
anova(gam1a)
newdat1a <- expand.grid(SI_LONG=agglarge$SI_LONG[1],SI_LATI=agglarge$SI_LATI[1], year=min(agglarge$year):max(agglarge$year), LE_GEAR_INNOV=agglarge$LE_GEAR_INNOV[1]) 
newdat1a$pred <- predict(gam1a, newdata=newdat1a, type="response" )
plot(exp(pred)~year, data=newdat1a)

gam1b       <- gam(value ~ te(SI_LONG,SI_LATI, year,k=5)+ LE_GEAR_INNOV,weights=TRIP, data=agglarge); 
gam.check(gam1b)
summary(gam1b)
plot(gam1b, all.terms=T)
anova(gam1b)
uniquerec <- agglarge[!duplicated(agglarge$LE_RECT),c("LE_RECT","SI_LONG","SI_LATI")]
newdat1b <- expand.grid(LE_RECT = uniquerec$LE_RECT, year=min(agglarge$year):max(agglarge$year), LE_GEAR_INNOV=agglarge$LE_GEAR_INNOV[1]) 
newdat1b <- merge(newdat1b,uniquerec)
newdat1b$pred <- predict(gam1b, newdata=newdat1b, type="response" )
plot(exp(pred)~year, data=newdat1b)
plot(exp(pred)~year, data=newdat1b[newdat1b$LE_RECT %in% lt4misdata, ])
meannewdat1b <- aggregate(pred~year, dat=newdat1b[newdat1b$LE_RECT %in% lt4misdata, ], FUN="mean")
plot(exp(pred)~year, data=meannewdat1b, ylim=c(0,0.12))

gam1c       <- gam(value ~ te(SI_LONG,SI_LATI,k=5)+te(year,k=10) + LE_GEAR_INNOV,weights=TRIP, data=agglarge); 
gam.check(gam1c)
summary(gam1c)
plot(gam1c, all.terms=T)
anova(gam1c)

gam1d       <- gam(value ~ te(SI_LONG,SI_LATI, by=as.factor(year),k=5)+ as.factor(year) + LE_GEAR_INNOV,weights=TRIP, data=agglarge); 
#gam.check(gam1d, type="")
gam.check(gam1d)
summary(gam1d)
plot(gam1d, all.terms=T)
anova(gam1d)
newdat1d <- expand.grid(LE_RECT = uniquerec$LE_RECT, year=min(agglarge$year):max(agglarge$year), LE_GEAR_INNOV=agglarge$LE_GEAR_INNOV[1]) 
newdat1d <- merge(newdat1d,uniquerec)
newdat1d$pred <- predict(gam1d, newdata=newdat1d, type="response" )
plot(exp(pred)~year, data=newdat1d)
plot(exp(pred)~year, data=newdat1d[newdat1d$LE_RECT %in% lt4misdata, ])

meannewdat1d <- aggregate(pred~year, dat=newdat1d[newdat1d$LE_RECT %in% lt4misdata, ], FUN="mean")
plot(exp(pred)~year, data=meannewdat1d, ylim=c(0,0.12))


AIC(gam1a, gam1b, gam1c, gam1d)
AIC(gam1a, gam1b, gam1c, gam1d, k= log(nobs(gam1a)))

origLpUE <- c(52.081,53.917,51.878,52.438,60.405,65.718,66.545,68.835,70.225,67.795,69.505,89.185,102.285,105.585,86.595,97.285,93.515,105.955,89.775,94.895,108.180)/1471

newdat1a <- expand.grid(SI_LONG=agglarge$SI_LONG[1],SI_LATI=agglarge$SI_LATI[1], year=min(agglarge$year):max(agglarge$year), LE_GEAR_INNOV=agglarge$LE_GEAR_INNOV[1])
newdat1a$pred <- predict(gam1a, newdata=newdat1a, type="response" )
plot(exp(pred)~year, data=newdat1a,xlab="Years",ylab="Predicted LpUE (KG/KWdays)",ylim=c(0,0.10),las=1,font.lab=2,type="b",pch="A")
grid();
lines(exp(pred)~year, data=newdat1a,type="b",col=1)
uniquerec <- agglarge[!duplicated(agglarge$LE_RECT),c("LE_RECT","SI_LONG","SI_LATI")]
newdat1b <- expand.grid(LE_RECT = uniquerec$LE_RECT, year=min(agglarge$year):max(agglarge$year), LE_GEAR_INNOV=agglarge$LE_GEAR_INNOV[1])
newdat1b <- merge(newdat1b,uniquerec)
newdat1b$pred <- predict(gam1b, newdata=newdat1b, type="response" )
meannewdat1b <- aggregate(pred~year, dat=newdat1b[newdat1b$LE_RECT %in% lt4misdata, ], FUN="mean")
lines(exp(pred)~year, data=meannewdat1b, col=2,type="b",pch="B")

newdat1c <- expand.grid(SI_LONG=agglarge$SI_LONG[1],SI_LATI=agglarge$SI_LATI[1], year=min(agglarge$year):max(agglarge$year), LE_GEAR_INNOV=agglarge$LE_GEAR_INNOV[1])
newdat1c$pred <- predict(gam1c, newdata=newdat1c, type="response" )
lines(exp(pred)~year, data=newdat1c,col=3,type="b",pch="C")

newdat1d <- expand.grid(LE_RECT = uniquerec$LE_RECT, year=min(agglarge$year):max(agglarge$year), LE_GEAR_INNOV=agglarge$LE_GEAR_INNOV[1])
newdat1d <- merge(newdat1d,uniquerec)
newdat1d$pred <- predict(gam1d, newdata=newdat1d, type="response" )
meannewdat1d <- aggregate(pred~year, dat=newdat1d[newdat1d$LE_RECT %in% lt4misdata, ], FUN="mean")
lines(exp(pred)~year, data=meannewdat1d, col=4,type="b",pch="D")

lines(y=origLpUE,x=1996:2016,col=5,pch="O",type="b")
legend("topleft",box.lty=0,col=1:5,pch=c("A","B","C","D","O"),legend=c(paste0("Model ",c("A","B","C","D")),"Original"),lty=1,lwd=2)

library(knitr)
tab <- cbind(newdat1a$year,exp(newdat1a$pred),exp(meannewdat1b$pred),exp(newdat1c$pred),exp(meannewdat1d$pred))
colnames(tab) <- c("Year",paste("Model",c("A","B","C","D")))
kable(tab)


for (var.plot in c("year","SI_LATI","SI_LONG")) {
#  png(file=paste(outPath,paste(var.plot),"wT.png",sep=""), bg="white",width = (25/2.54), height = (25/2.54), units = "in", pointsize = 12,res=100)
#  par(mfrow=c(2,2),mar=c(2, 2, 2, 1),oma=c(0,0,0,0))
  for (gr in c(1:length(grs))) {
    
    agglarge.sel<-agglarge[agglarge$LE_GEAR_INNOV==grs[gr],]
    print(grs[gr])
    print(exp(weighted.mean(agglarge.sel$value, agglarge.sel$TRIP)))
    #gam1       <- gam(value ~ as.factor(year)+s(SI_LONG,k=4)+s(SI_LATI,k=4),weights=TRIP, data=agglarge.sel); anova(gam1)
    gam1       <- gam(value ~ s(SI_LONG,k=4)+s(year,k=4)+s(SI_LATI,k=4),weights=TRIP, data=agglarge.sel); anova(gam1)
    vars.lab<-data.frame(var.plot=c("year",
                                    "SI_LATI",
                                    "SI_LONG"),
                         
                         xlab.title=c("",
                                      "",
                                      ""))
    pred.dat<-as.data.frame(matrix(NA,nrow=100,ncol=nrow(vars.lab)))
    names(pred.dat)<-as.character(vars.lab$var.plot)
    #var.plot<-"year"; 
    if (var.plot == "year")  
    {pred.dat$year<-rep(unique(agglarge.sel$year),20)[1:100]
    for (K in 2:nrow(vars.lab))
    {
      pred.dat[,K]<-rep(mean(agglarge.sel[,as.character(vars.lab$var.plot[K])],na.rm=T),100)
    }
    }
    if (var.plot == "SI_LATI")  
    {pred.dat$SI_LATI<-rep(unique(agglarge.sel$SI_LATI),20)[1:100]
    pred.dat[,1]<-2016
    pred.dat[,3]<-rep(mean(agglarge.sel[,as.character(vars.lab$var.plot[3])],na.rm=T),100)
    pred.dat[,3]<-4}
    if (var.plot == "SI_LONG")  
    {pred.dat$SI_LONG<-rep(unique(agglarge.sel$SI_LONG),20)[1:100]
    pred.dat[,1]<-2016
    pred.dat[,2]<-rep(mean(agglarge.sel[,as.character(vars.lab$var.plot[2])],na.rm=T),100)
    }
    
    head(pred.dat)
    pred.dat[,var.plot]<-seq(min(agglarge.sel[,var.plot]),max(agglarge[,var.plot]),length=100)
    head(pred.dat)
    # Make plot
    f<-function(x){exp(x)}
    prediction<-predict(gam1,newdata=pred.dat,se=T)
    head(prediction)
    ylims=c(0,max(exp(prediction$fit+prediction$se.fit)))
    #ylabs="year"
    ylims=c(0,0.20)
    if (var.plot == "year")  ylims=c(0,0.120)
    pred.dat$prediction<-prediction$fit
    pred.dat$se<-prediction$se.fit
    plot(pred.dat[,var.plot],f(pred.dat$prediction),type="l",col="red",main=paste(grs[gr]),
         lwd=2,ylim=ylims,xlab=var.plot,ylab="LpUE")
    rug(agglarge[,var.plot])
    lines(pred.dat[,var.plot],f(pred.dat$prediction+pred.dat$se*2),lty=2)
    lines(pred.dat[,var.plot],f(pred.dat$prediction-pred.dat$se*2),lty=2)
  }
#  dev.off()
}




