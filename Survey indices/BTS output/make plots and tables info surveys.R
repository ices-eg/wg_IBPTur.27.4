library(sas7bdat); library(FLCore); library(grid)

snsagesamp <- read.sas7bdat( "w://IMARES/Data/ICES-WG/IBPTURBOT/2017/survey indices/SNS sas files/age_tur.sas7bdat")

table(snsagesamp[snsagesamp$aged=="aged" & snsagesamp$PGM_CODE=="SNS",]$year, snsagesamp[snsagesamp$aged=="aged" & snsagesamp$PGM_CODE=="SNS" ,]$year-snsagesamp[snsagesamp$aged=="aged" &snsagesamp$PGM_CODE == "SNS" ,]$yc)
table(snsagesamp[snsagesamp$aged=="aged"& snsagesamp$PGM_CODE=="SNS",]$year)

snsnjlen <- read.sas7bdat( "w://IMARES/Data/ICES-WG/IBPTURBOT/2017/survey indices/SNS sas files/snsnjlen_tur.sas7bdat")

snsnjlen[is.nan(snsnjlen$snumtot),]$snumtot <- 0
table(snsnjlen$year,snsnjlen$snumtot)

aggregate(snumtot~ year, data= snsnjlen, FUN= "sum")

table(snsnjlen$year)

###############################################
# BTS info
###############################################
btsagesamp <- read.sas7bdat( "w://IMARES/Data/ICES-WG/IBPTURBOT/2017/survey indices/BTS sas files/age_tur.sas7bdat")

table(btsagesamp[btsagesamp$aged=="aged" & btsagesamp$PGM_CODE=="BTS",]$year, btsagesamp[btsagesamp$aged=="aged" & btsagesamp$PGM_CODE=="BTS" ,]$year-btsagesamp[btsagesamp$aged=="aged" &btsagesamp$PGM_CODE == "BTS" ,]$yc)
table(btsagesamp[btsagesamp$aged=="aged" & btsagesamp$PGM_CODE=="BTS",]$year)

btsnjlen <- read.sas7bdat( "w://IMARES/Data/ICES-WG/IBPTURBOT/2017/survey indices/BTS sas files/btsnjlen_tur.sas7bdat")

btsnjlen[is.nan(btsnjlen$snumtot),]$snumtot <- 0
btsnjlen[btsnjlen$year %in% 1986:1997, ]$snumtot <- btsnjlen[btsnjlen$year %in% 1986:1997, ]$snumtot * 2

table(btsnjlen$year,btsnjlen$snumtot)
tail(table(btsnjlen$year,btsnjlen$snumtot),6)

table(btsnjlen[btsnjlen$snumtot>=15,]$year)

aggregate(snumtot~ year, data= btsnjlen, FUN= "sum")

table(btsnjlen$year)


###########################################################
# READ LOWESTOFT FILES AND PLOT consistencies SNS and BTS-ISIS for different starting years
############################################################

lowestoftPath <- "w://IMARES/Data/ICES-WG/IBPTURBOT/2017/assessment runs/lowestoft files/"

indices            <- readFLIndices(paste(lowestoftPath,"fleet_tmb_longer_LPUE.txt", sep=""), na.strings="-1")

names(indices)

indices <- FLIndices(list(indices[[1]],indices[[2]]))

indexVals <- lapply(indices, index)
numIndices <- length(indexVals)
indMPs <- list()

plot(indices[[1]])

#### Run indices diagnostics
indices.short <- FLIndices(list(window(indices[[1]],2004, 2016),window(indices[[2]],1991,2016)))

idxcrop <-  indices
indsN01 <- FLQuants(lapply( mcf( lapply(idxcrop, index)), function(x){x <- FLQuant(aperm(apply(x@.Data, c(1,3,4,5,6), scale),c(2,1,3,4,5,6)), dimnames= dimnames(x))}))

names(indsN01)   <- names(indices)
akey             <- simpleKey(text=names(indsN01), points=F, lines=F, columns=2, cex=1.5, col=c("red","black","blue","gray","orange","magenta"))
#akey$text$lab[1] <- "BTS-ISIS"; akey$text$lab[2] <- "SNS"
xyplot(data~year | factor(age), data=indsN01, type="b", key=akey, groups=qname, pch=19, 
       col=c("red","black","blue","gray","orange","magenta"),as.table=TRUE, scales="free",  layout=c(3,3), xlim=c(1985,2016), ylim=c(-1.5,3.5))

idxcrop <-  indices.short
indsN01 <- FLQuants(lapply( mcf( lapply(idxcrop, index)), function(x){x <- FLQuant(aperm(apply(x@.Data, c(1,3,4,5,6), scale),c(2,1,3,4,5,6)), dimnames= dimnames(x))}))

names(indsN01)   <- names(indices.short)
akey             <- simpleKey(text=names(indsN01), points=F, lines=F, columns=2, cex=1.5, col=c("red","black","blue","gray","orange","magenta"))
#akey$text$lab[1] <- "BTS-ISIS"; akey$text$lab[2] <- "SNS"
xyplot(data~year | factor(age), data=indsN01, type="b", key=akey, groups=qname, pch=19, 
       col=c("red","black","blue","gray","orange","magenta"),as.table=TRUE, scales="free",  layout=c(3,3), xlim=c(1985,2016), ylim=c(-1.5,3.5))


for (i in c(1,2)) {
  plot(indices[[i]], type="internal", main=paste0(name(indices[[i]])," ", range(indices[[i]])[["minyear"]],"-",range(indices[[i]])[["maxyear"]]))
}


for (i in c(1,2)) {
  plot(indices.short[[i]], type="internal", main=paste0(name(indices.short[[i]])," ", range(indices.short[[i]])[["minyear"]],"-",range(indices.short[[1]])[["maxyear"]]))
}






