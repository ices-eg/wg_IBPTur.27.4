# Comparison of IBP 2017 with WGNSSK and IBP2018 runs
rm(list=ls())

# LOAD FILES
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/baserun/2017_final_finalassessmentOut.RData")
T1 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/baserun/TUR.sam_2018_assesment_SSB.Rdata")
T2 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/baserun/TUR.sam_2018_assesment_EB.Rdata")
T3 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/baserun/TUR.sam_2018_base_run_SSB.Rdata")
T4  <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/baserun/TUR.sam_2018_base_run_EB.Rdata")
T5 <- TUR.sam

plot(ssb(T1)$value,x=ssb(T1)$year,col=3,lwd=4,type="l",ylim=c(0,18000),las=1, main= "SSB",ylab="SSB (tonnes)",xlab=" ")
lines(y=ssb(T2)$value, x=ssb(T2)$year ,type="l",lwd=2)
lines(y=ssb(T3)$value, x=ssb(T3)$year ,type="l",las=1,lwd=2,col=2)
lines(y=ssb(T4)$value, x=ssb(T4)$year ,type="l",las=1,lwd=2,col=4)
lines(y=ssb(T5)$value, x=ssb(T5)$year ,type="l",las=1,lwd=2,col=5)
legend(1982, 5000, legend=c("Final_IBP_2017_SSB", "Final_2018_SSB","Final_2018_EB", "Base_2018_SSB","Base_2018_EB"),
       col=c("3", "1","2","4","5"), lty=1, cex=0.8)

plot(fbar(T1)$value,x=fbar(T1)$year,col=3,lwd=4,type="l",ylim=c(0,1),las=1, main= "fbar",ylab="fbar",xlab=" ")
lines(y=fbar(T2)$value, x=fbar(T2)$year ,type="l",lwd=2)
lines(y=fbar(T3)$value, x=fbar(T3)$year ,type="l",las=1,lwd=2,col=2)
lines(y=fbar(T4)$value, x=fbar(T4)$year ,type="l",las=1,lwd=2,col=4)
lines(y=fbar(T5)$value, x=fbar(T5)$year ,type="l",las=1,lwd=2,col=5)
legend(1982, 0.2, legend=c("Final_IBP_2017_SSB", "Final_2018_SSB","Final_2018_EB", "Base_2018_SSB","Base_2018_EB"),
       col=c("3", "1","2","4","5"), lty=1, cex=0.8)

plot(rec(T1)$value,x=rec(T1)$year,col=3,lwd=4,type="l",ylim=c(0,10000),las=1, main= "rec",ylab="rec (tonnes)",xlab=" ")
lines(y=rec(T2)$value, x=rec(T2)$year ,type="l",lwd=2)
lines(y=rec(T3)$value, x=rec(T3)$year ,type="l",las=1,lwd=2,col=2)
lines(y=rec(T4)$value, x=rec(T4)$year ,type="l",las=1,lwd=2,col=4)
lines(y=rec(T5)$value, x=rec(T5)$year ,type="l",las=1,lwd=2,col=5)
legend(1982, 1800, legend=c("Final_IBP_2017_SSB", "Final_2018_SSB","Final_2018_EB", "Base_2018_SSB","Base_2018_EB"),
       col=c("3", "1","2","4","5"), lty=1, cex=0.8)



# MAKE PLOT WITH ALL SENSITIVITY RUNS TOGETHER
# LOAD FILES

# F.states
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/01_Ctrl.States/Sens_Ctrl.States_1assessmentOut.RData")
ST1 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/01_Ctrl.States/Sens_Ctrl.States_2assessmentOut.RData")
ST2 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/01_Ctrl.States/Sens_Ctrl.States_3assessmentOut.RData")
ST3 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/01_Ctrl.States/Sens_Ctrl.States_4assessmentOut.RData")
ST4  <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/01_Ctrl.States/Sens_Ctrl.States_5assessmentOut.RData")
ST5 <- TUR.sam

plot(ssb(ST1)$value,x=ssb(ST1)$year,col=3,lwd=4,type="l",ylim=c(0,18000),las=1, main= "SSB",ylab="SSB (tonnes)",xlab=" ")
lines(y=ssb(ST2)$value, x=ssb(ST2)$year ,type="l",lwd=2)
lines(y=ssb(ST3)$value, x=ssb(ST3)$year ,type="l",las=1,lwd=2,col=2)
lines(y=ssb(ST4)$value, x=ssb(ST4)$year ,type="l",las=1,lwd=2,col=4)
lines(y=ssb(ST5)$value, x=ssb(ST5)$year ,type="l",las=1,lwd=2,col=5)
legend(1982, 5000, legend=c("State_1", "State_2","State_3", "State_4","State_5"),
       col=c("3", "1","2","4","5"), lty=1, cex=0.8)

# Cor.F
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/02_Cor.F/Sens_Cor.F_0assessmentOut.RData")
CF0 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/02_Cor.F/Sens_Cor.F_1assessmentOut.RData")
CF1 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/02_Cor.F/Sens_Cor.F_2assessmentOut.RData")
CF2 <- TUR.sam

plot(ssb(CF0)$value,x=ssb(CF0)$year,col=3,lwd=4,type="l",ylim=c(0,18000),las=1, main= "SSB",ylab="SSB (tonnes)",xlab=" ")
lines(y=ssb(CF1)$value, x=ssb(CF1)$year ,type="l",lwd=2)
lines(y=ssb(CF2)$value, x=ssb(CF2)$year ,type="l",las=1,lwd=2,col=2)
legend(1982, 5000, legend=c("Cor.F_0", "Cor.F_1","Cor.F_2"),
       col=c("3", "1","2"), lty=1, cex=0.8)

# SNS.bind
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/03_SNS.bind/Sens_SNS.bind_1assessmentOut.RData")
SNS1 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/03_SNS.bind/Sens_SNS.bind_2assessmentOut.RData")
SNS2 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/03_SNS.bind/Sens_SNS.bind_3assessmentOut.RData")
SNS3 <- TUR.sam

plot(ssb(SNS1)$value,x=ssb(SNS1)$year,col=3,lwd=4,type="l",ylim=c(0,18000),las=1, main= "SSB",ylab="SSB (tonnes)",xlab=" ")
lines(y=ssb(SNS2)$value, x=ssb(SNS2)$year ,type="l",lwd=2)
lines(y=ssb(SNS3)$value, x=ssb(SNS3)$year ,type="l",las=1,lwd=2,col=2)
legend(1982, 5000, legend=c("SNS.bind1", "SNS.bind2","SNS.bind3"),
       col=c("3", "1","2"), lty=1, cex=0.8)

# BTS.bind
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/04_BTS.bind/Sens_BTS.bind_1assessmentOut.RData")
BTS1 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/04_BTS.bind/Sens_BTS.bind_2assessmentOut.RData")
BTS2 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/04_BTS.bind/Sens_BTS.bind_3assessmentOut.RData")
BTS3 <- TUR.sam

plot(ssb(BTS1)$value,x=ssb(BTS1)$year,col=3,lwd=4,type="l",ylim=c(0,18000),las=1, main= "SSB",ylab="SSB (tonnes)",xlab=" ")
lines(y=ssb(BTS2)$value, x=ssb(BTS2)$year ,type="l",lwd=2)
lines(y=ssb(BTS3)$value, x=ssb(BTS3)$year ,type="l",las=1,lwd=2,col=2)
legend(1982, 5000, legend=c("BTS.bind1", "BTS.bind2","BTS.bind3"),
       col=c("3", "1","2"), lty=1, cex=0.8)

#f.vars
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/05_f.vars/sens_f.vars_1assessmentOut.RData")
Fvar1 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/05_f.vars/sens_f.vars_2assessmentOut.RData")
Fvar2 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/05_f.vars/sens_f.vars_4assessmentOut.RData")
Fvar4 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/05_f.vars/sens_f.vars_6assessmentOut.RData")
Fvar6 <- TUR.sam

# obs.vars Catch
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/06_Obs.vars._Catch/Sens_Obs.vars_catch_1assessmentOut.RData")
OVCat1 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/06_Obs.vars._Catch/Sens_Obs.vars_catch_2assessmentOut.RData")
OVCat2 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/06_Obs.vars._Catch/Sens_Obs.vars_catch_3assessmentOut.RData")
OVCat3 <- TUR.sam

# obs.vars SNS
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/07_Obs.vars_SNS/Sens_Obs.vars_SNS_1assessmentOut.RData")
OVSNS1 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/07_Obs.vars_SNS/Sens_Obs.vars_SNS_2assessmentOut.RData")
OVSNS2 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/07_Obs.vars_SNS/Sens_Obs.vars_SNS_3assessmentOut.RData")
OVSNS3 <- TUR.sam

# obs.vars BTS
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/08_Obs.vars_BTS/Sens_Obs.vars_BTS_1assessmentOut.RData")
OVBTS1 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/08_Obs.vars_BTS/Sens_Obs.vars_BTS_2assessmentOut.RData")
OVBTS2 <- TUR.sam
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/08_Obs.vars_BTS/Sens_Obs.vars_BTS_3assessmentOut.RData")
OVBTS3 <- TUR.sam

#cor.obs SNS
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/09 Cor.Obs_SNS/Sens_Cor.Obs_SNS_1assessmentOut.RData")
COSNS <- TUR.sam

# final
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/10_Final/IBP_final_run__1_assessmentOut.RData")
Final <- TUR.sam

outPath <- "D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/plots/"

# SSB
png(file= paste0(outPath,"All_runs_SSB.png"),units="px",pointsize = 12,width=960,height = 600)
plot(ssb(ST1)$value,x=ssb(ST1)$year,col=1,lwd=1,type="l",ylim=c(0,16000),las=1, ylab="SSB (tonnes)",
     xlab=" ",yaxs="i",cex.lab=1.3)
lines(y=ssb(ST2)$value, x=ssb(ST2)$year ,lwd=1,col=1)
lines(y=ssb(ST3)$value, x=ssb(ST3)$year ,lwd=1,col=1)
lines(y=ssb(ST4)$value, x=ssb(ST4)$year ,lwd=1,col=1)
lines(y=ssb(ST5)$value, x=ssb(ST5)$year ,lwd=1,col=1)
lines(y=ssb(CF0)$value, x=ssb(CF0)$year ,lwd=1,col="grey")
lines(y=ssb(CF1)$value, x=ssb(CF1)$year ,lwd=1,col="grey")
lines(y=ssb(CF2)$value, x=ssb(CF2)$year ,lwd=1,col="grey")
lines(y=ssb(SNS1)$value, x=ssb(SNS1)$year ,lwd=1,col="grey15")
lines(y=ssb(SNS2)$value, x=ssb(SNS2)$year ,lwd=1,col="grey15")
lines(y=ssb(SNS3)$value, x=ssb(SNS3)$year ,lwd=1,col="grey15")
lines(y=ssb(BTS1)$value, x=ssb(BTS1)$year ,lwd=1,col="grey35")
lines(y=ssb(BTS2)$value, x=ssb(BTS2)$year ,lwd=1,col="grey35")
lines(y=ssb(BTS3)$value, x=ssb(BTS3)$year ,lwd=1,col="grey35")
lines(y=ssb(Fvar1)$value, x=ssb(Fvar1)$year ,lwd=1,col="grey55")
lines(y=ssb(Fvar2)$value, x=ssb(Fvar2)$year ,lwd=1,col="grey55")
lines(y=ssb(Fvar4)$value, x=ssb(Fvar4)$year ,lwd=1,col="grey55")
lines(y=ssb(Fvar6)$value, x=ssb(Fvar6)$year ,lwd=1,col="grey55")
lines(y=ssb(OVCat1)$value, x=ssb(OVCat1)$year ,lwd=1,col="grey75")
lines(y=ssb(OVCat2)$value, x=ssb(OVCat2)$year ,lwd=1,col="grey75")
lines(y=ssb(OVCat3)$value, x=ssb(OVCat3)$year ,lwd=1,col="grey75")
lines(y=ssb(OVSNS1)$value, x=ssb(OVSNS1)$year ,lwd=1,col="grey95")
lines(y=ssb(OVSNS2)$value, x=ssb(OVSNS2)$year ,lwd=1,col="grey95")
lines(y=ssb(OVSNS3)$value, x=ssb(OVSNS3)$year ,lwd=1,col="grey95")
lines(y=ssb(OVBTS1)$value, x=ssb(OVBTS1)$year ,lwd=1,col="grey5")
lines(y=ssb(OVBTS2)$value, x=ssb(OVBTS2)$year ,lwd=1,col="grey5")
lines(y=ssb(OVBTS3)$value, x=ssb(OVBTS3)$year ,lwd=1,col="grey5")
lines(y=ssb(COSNS)$value, x=ssb(COSNS)$year ,lwd=1,col="grey45")
lines(y=ssb(Final)$value, x=ssb(Final)$year ,lwd=1,col="red")
abline(h=10000,col=2,lty=2)
dev.off()
# fbar

png(file= paste0(outPath,"All_runs_fbar.png"),units="px",pointsize = 12,width=960,height = 600)
plot(fbar(ST1)$value,x=fbar(ST1)$year,col=1,lwd=1,type="l",ylim=c(0,1),las=1, ylab="fbar (ages 2-6)",
     xlab=" ",yaxs="i",cex.lab=1.3)
lines(y=fbar(ST2)$value, x=fbar(ST2)$year ,lwd=1,col=1)
lines(y=fbar(ST3)$value, x=fbar(ST3)$year ,lwd=1,col=1)
lines(y=fbar(ST4)$value, x=fbar(ST4)$year ,lwd=1,col=1)
lines(y=fbar(ST5)$value, x=fbar(ST5)$year ,lwd=1,col=1)
lines(y=fbar(CF0)$value, x=fbar(CF0)$year ,lwd=1,col="grey")
lines(y=fbar(CF1)$value, x=fbar(CF1)$year ,lwd=1,col="grey")
lines(y=fbar(CF2)$value, x=fbar(CF2)$year ,lwd=1,col="grey")
lines(y=fbar(SNS1)$value, x=fbar(SNS1)$year ,lwd=1,col="grey15")
lines(y=fbar(SNS2)$value, x=fbar(SNS2)$year ,lwd=1,col="grey15")
lines(y=fbar(SNS3)$value, x=fbar(SNS3)$year ,lwd=1,col="grey15")
lines(y=fbar(BTS1)$value, x=fbar(BTS1)$year ,lwd=1,col="grey35")
lines(y=fbar(BTS2)$value, x=fbar(BTS2)$year ,lwd=1,col="grey35")
lines(y=fbar(BTS3)$value, x=fbar(BTS3)$year ,lwd=1,col="grey35")
lines(y=fbar(Fvar1)$value, x=fbar(Fvar1)$year ,lwd=1,col="grey55")
lines(y=fbar(Fvar2)$value, x=fbar(Fvar2)$year ,lwd=1,col="grey55")
lines(y=fbar(Fvar4)$value, x=fbar(Fvar4)$year ,lwd=1,col="grey55")
lines(y=fbar(Fvar6)$value, x=fbar(Fvar6)$year ,lwd=1,col="grey55")
lines(y=fbar(OVCat1)$value, x=fbar(OVCat1)$year ,lwd=1,col="grey75")
lines(y=fbar(OVCat2)$value, x=fbar(OVCat2)$year ,lwd=1,col="grey75")
lines(y=fbar(OVCat3)$value, x=fbar(OVCat3)$year ,lwd=1,col="grey75")
lines(y=fbar(OVSNS1)$value, x=fbar(OVSNS1)$year ,lwd=1,col="grey95")
lines(y=fbar(OVSNS2)$value, x=fbar(OVSNS2)$year ,lwd=1,col="grey95")
lines(y=fbar(OVSNS3)$value, x=fbar(OVSNS3)$year ,lwd=1,col="grey95")
lines(y=fbar(OVBTS1)$value, x=fbar(OVBTS1)$year ,lwd=1,col="grey5")
lines(y=fbar(OVBTS2)$value, x=fbar(OVBTS2)$year ,lwd=1,col="grey5")
lines(y=fbar(OVBTS3)$value, x=fbar(OVBTS3)$year ,lwd=1,col="grey5")
lines(y=fbar(COSNS)$value, x=fbar(COSNS)$year ,lwd=1,col="grey45")
lines(y=fbar(Final)$value, x=fbar(Final)$year ,lwd=1,col="red")
#abline(h=10000,col=2,lty=2)
dev.off()

# rec
png(file= paste0(outPath,"All_runs_Rec.png"),units="px",pointsize = 12,width=960,height = 600)
plot(rec(ST1)$value,x=rec(ST1)$year,col=1,lwd=1,type="l",ylim=c(0,10000),las=1, 
     ylab="recruitement (thousands)",xlab=" ",yaxs="i",cex.lab=1.3)
lines(y=rec(ST2)$value, x=rec(ST2)$year ,lwd=1,col=1)
lines(y=rec(ST3)$value, x=rec(ST3)$year ,lwd=1,col=1)
lines(y=rec(ST4)$value, x=rec(ST4)$year ,lwd=1,col=1)
lines(y=rec(ST5)$value, x=rec(ST5)$year ,lwd=1,col=1)
lines(y=rec(CF0)$value, x=rec(CF0)$year ,lwd=1,col="grey")
lines(y=rec(CF1)$value, x=rec(CF1)$year ,lwd=1,col="grey")
lines(y=rec(CF2)$value, x=rec(CF2)$year ,lwd=1,col="grey")
lines(y=rec(SNS1)$value, x=rec(SNS1)$year ,lwd=1,col="grey15")
lines(y=rec(SNS2)$value, x=rec(SNS2)$year ,lwd=1,col="grey15")
lines(y=rec(SNS3)$value, x=rec(SNS3)$year ,lwd=1,col="grey15")
lines(y=rec(BTS1)$value, x=rec(BTS1)$year ,lwd=1,col="grey35")
lines(y=rec(BTS2)$value, x=rec(BTS2)$year ,lwd=1,col="grey35")
lines(y=rec(BTS3)$value, x=rec(BTS3)$year ,lwd=1,col="grey35")
lines(y=rec(Fvar1)$value, x=rec(Fvar1)$year ,lwd=1,col="grey55")
lines(y=rec(Fvar2)$value, x=rec(Fvar2)$year ,lwd=1,col="grey55")
lines(y=rec(Fvar4)$value, x=rec(Fvar4)$year ,lwd=1,col="grey55")
lines(y=rec(Fvar6)$value, x=rec(Fvar6)$year ,lwd=1,col="grey55")
lines(y=rec(OVCat1)$value, x=rec(OVCat1)$year ,lwd=1,col="grey75")
lines(y=rec(OVCat2)$value, x=rec(OVCat2)$year ,lwd=1,col="grey75")
lines(y=rec(OVCat3)$value, x=rec(OVCat3)$year ,lwd=1,col="grey75")
lines(y=rec(OVSNS1)$value, x=rec(OVSNS1)$year ,lwd=1,col="grey95")
lines(y=rec(OVSNS2)$value, x=rec(OVSNS2)$year ,lwd=1,col="grey95")
lines(y=rec(OVSNS3)$value, x=rec(OVSNS3)$year ,lwd=1,col="grey95")
lines(y=rec(OVBTS1)$value, x=rec(OVBTS1)$year ,lwd=1,col="grey5")
lines(y=rec(OVBTS2)$value, x=rec(OVBTS2)$year ,lwd=1,col="grey5")
lines(y=rec(OVBTS3)$value, x=rec(OVBTS3)$year ,lwd=1,col="grey5")
lines(y=rec(COSNS)$value, x=rec(COSNS)$year ,lwd=1,col="grey45")
lines(y=rec(Final)$value, x=rec(Final)$year ,lwd=1,col="red")
dev.off()


# Comparison between final run IBP 2017 and IBP 2018
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/baserun/2017_final_finalassessmentOut.RData")
IBP2017 <- TUR.sam
# final
load("D:/Repository/wg_IBPTur.27.4/IBP_2018/Trial_runs_2018/Output/parSettings/10_Final/IBP_final_run__1_assessmentOut.RData")
IBP2018 <- TUR.sam

png(file= paste0(outPath,"Comp_fin_runs_IBP_SSB.png"),units="px",pointsize = 12,width=960,height = 600)
plot(ssb(IBP2017)$value,x=ssb(IBP2017)$year,col="gray25",lwd=2,type="l",ylim=c(0,20000),las=1, ylab="ssb (tonnes)",
     xlab=" ",yaxs="i",cex.lab=1.3)
lines(y=ssb(IBP2018)$value, x=ssb(IBP2018)$year ,lwd=2,col="gray55")
legend(2010, 18000, legend=c("IBP_2017", "IBP_2018"),
       col=c("gray25","gray55"), lty=1, cex=1.2,lwd=2)
dev.off()

png(file= paste0(outPath,"Comp_fin_runs_IBP_fbar.png"),units="px",pointsize = 12,width=960,height = 600)
plot(fbar(IBP2017)$value,x=fbar(IBP2017)$year,col="gray25",lwd=2,type="l",ylim=c(0,1),las=1, ylab="fbar (age 2-6)",
     xlab=" ",yaxs="i",cex.lab=1.3)
lines(y=fbar(IBP2018)$value, x=fbar(IBP2018)$year ,lwd=2,col="gray55")
legend(2010, 0.9, legend=c("IBP_2017", "IBP_2018"),
       col=c("gray25","gray55"), lty=1, cex=1.2,lwd=2)
dev.off()

png(file= paste0(outPath,"Comp_fin_runs_IBP_rec.png"),units="px",pointsize = 12,width=960,height = 600)
plot(rec(IBP2017)$value,x=rec(IBP2017)$year,col="gray25",lwd=2,type="l",ylim=c(0,10000),las=1, ylab="recruitement (thousands)",
     xlab=" ",yaxs="i",cex.lab=1.3)
lines(y=rec(IBP2018)$value, x=rec(IBP2018)$year ,lwd=2,col="gray55")
legend(1982, 9000, legend=c("IBP_2017", "IBP_2018"),
       col=c("gray25","gray55"), lty=1, cex=1.2,lwd=2)
dev.off()