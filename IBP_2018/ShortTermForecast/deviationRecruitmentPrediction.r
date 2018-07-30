store <- matrix(NA,nrow=length(2003:2017),ncol=1,dimnames=list(year=2003:2017,1))
for(iYr in 2003:2017)
  store[ac(iYr),1] <- mean(rec(TUR[,ac(iYr:2017)]))
  
  
plot(y=store[,1],x=15:1,xlab="Geomean period length",ylab="Recruitment",type="b")

rec <- rec(TUR)[,drop=T]

storeDev <- matrix(NA,nrow=length(2006:2016),ncol=6,dimnames=list(AssYear=2006:2016,geoyears=c(3,5,10,15,20,"entire")))
for(iYr in 2006:2016){
  year3   <- exp(mean(log(rec[ac((iYr-2):iYr)])))
  year5   <- exp(mean(log(rec[ac((iYr-4):iYr)])))
  year10  <- exp(mean(log(rec[ac((iYr-9):iYr)])))
  year15  <- exp(mean(log(rec[ac((iYr-14):iYr)])))
  year20  <- exp(mean(log(rec[ac((iYr-19):iYr)])))
  entire  <- exp(mean(log(rec)))


  pred    <- rec[ac(iYr+1)]
  
  dev3    <- (year3-pred)/pred
  dev5    <- (year5-pred)/pred
  dev10   <- (year10-pred)/pred
  dev15   <- (year15-pred)/pred
  dev20   <- (year20-pred)/pred
  entire  <- (entire-pred)/pred

  storeDev[ac(iYr),] <- c(dev3,dev5,dev10,dev15,dev20,entire)

}

matplot(y=storeDev*100,x=2006:2016,xlab="Assessment year",ylab="Deviation in %",pch=c("3","5","0","1","2","e"),type="b")
abline(h=0,col=4,lty=3)