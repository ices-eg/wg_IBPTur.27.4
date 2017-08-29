save(TUR.sam,TUR,TUR.tun,TUR.ctrl,TUR.retro,file=file.path(outPath,paste0(run,"_",sens,"assessmentOut.RData")))
source(paste0(codePath,"retroResidual.r"))

pdf(file.path(outPath,paste0(run,"_",sens,"assessmentOut.pdf")))
#png(file.path(outPath,"figures - %02d.png"),units = "px", height=800,width=672, bg = "white")
try(residual.diagnostics(TUR.sam))
print(cor.plot(TUR.sam))
obscv.plot(TUR.sam)
# figure - catchabilities at age from HERAS
catch <- catchabilities(TUR.sam)
print(xyplot(value+ubnd+lbnd ~ age | fleet,catch,
             scale=list(alternating=FALSE,y=list(relation="free")),as.table=TRUE,
             type="l",lwd=c(2,1,1),col=c("black","grey","grey"),
             subset=fleet %in% c("SNS","BTS-ISIS","NL_LPUE_age","BTS_DG"),
             main="Survey catchability parameters",ylab="Catchability",xlab="Age"))
obsvar.plot(TUR.sam)
# figure - fishing age selectivity per year
sel.pat <- merge(f(TUR.sam),fbar(TUR.sam),
               by="year",suffixes=c(".f",".fbar"))
sel.pat$sel <- sel.pat$value.f/sel.pat$value.fbar
sel.pat$age <- as.numeric(as.character(sel.pat$age))
print(xyplot(sel ~ age|sprintf("%i's",floor((year+2)/5)*5),sel.pat,
             groups=year,type="l",as.table=TRUE,
             scale=list(alternating=FALSE),
             main="Selectivity of the Fishery by Pentad",xlab="Age",ylab="F/Fbar"))
print(plot(TUR.sam))
dat <- subset(residuals(TUR.sam),fleet=="catch")
print(xyplot(age ~ year,data=dat,cex=dat$std.res,col="black",main="Residuals by year Catch",
panel=function(...){
    lst <- list(...)
    panel.xyplot(lst$x,lst$y,pch=ifelse(lst$cex>0,1,19),col="black",cex=abs(lst$cex))
}))

# figure - acosutic index residuals per year per age
dat <- subset(residuals(TUR.sam),fleet %in% c("SNS","BTS-ISIS","NL_LPUE_age"))
print(xyplot(age ~ year|fleet,data=dat,cex=dat$std.res,col="black",main="Residuals by survey",
panel=function(...){
    lst <- list(...)
    panel.xyplot(lst$x,lst$y,pch=ifelse(lst$cex>0,1,19),col="black",cex=abs(lst$cex))
}))

obscor.plot(TUR.sam)

print(plot(TUR.retro))
print(retroParams(TUR.retro))
print(retroSelectivity(TUR.retro,2007:2016))

par(mfrow=c(2,1))
print(plot(AIC(TUR.sam),ylab="AIC",xaxt="n",las=2,pch=19,xlab=""))
print(grid())
print(plot(mean(mohns.rho(TUR.retro,2016,7)[1:7,"rho"]),xlab="",ylab="Mohns rho (7-year peel)",xaxt="n",las=2,pch=19))
grid()
print(grid())

dev.off()