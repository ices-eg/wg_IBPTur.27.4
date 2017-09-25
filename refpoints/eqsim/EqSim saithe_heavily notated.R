### ---------------------------------------------
## setting constants
Fcv<-0.212
Fphi<-0.423
Blim<-113000     # Bloss
## need to force FIT to use Bloss, not estimated Blim from segreg -- use segreg3 where apply breakpt directly
Bpa<-round(Blim*1.4,-3) #exp(1.645*0.1524100),-3)       ##  in samtmb.RData sdrep$sd[96] #check sdrep$value[96] first 
Flim<- 0.712                                         #F50 from landings (not catch)
Fpa <- round(Flim/1.4,3) #exp(1.645 * 2.099004e-01),2)         # sdrep$sd[192]  = -1.517777
#Btrigger = 277000   # this is from the 5th percentil on BFMSY this is greater than Bpa & the current Btrigger
## estimating 5th percentile of SSB2014 
Btrigger<- round(round(exp(12.57701-(1.645*0.1524100))),-3)  ## 226000 ## this is the sd of SSB2014 sdrep$sd[96]


#############--- Simulations Full time series ---###################
###### no Fcv and Fphi ------------------ 
## (to get Flim--F that gives 50% probability of SSB>Blim)  -- use the F associated with the catch 
SIM101 <- eqsim_run(FIT,  bio.years = c(2005, 2014), bio.const = FALSE,sel.years = c(2005, 2014), sel.const = FALSE,Fcv=0, Fphi=0,Btrigger = 0,Blim=Blim,Bpa=Bpa,Fscan = seq(0,1.2,len=61),verbose=FALSE)

setwd("C:\\Users\\jenniferd\\Documents\\_2016\\WGNSSK\\WKNSEA\\Assessements\\final\\Eqsim\\")
png("plot_saithe_full time period no Btrigger no Fcv no Fphi.png")
eqsim_plot(SIM101,catch="FALSE")
dev.off()

Coby.fit(SIM101,outfile='saithe full time series no Btrigger Fcv=0 and Fphi=0')


##----------------------------------------
## gives FMSY -- should give 0.345 for median FMSY with new Blim and with discards included for age 4+ ------------------
## look at the F associated with the landings only
#range without advice rule (no Btrigger, etc) 0.198 lower, 0.574 upper
## Fp05 also from this sim
SIM1 <- eqsim_run(FIT,  bio.years = c(2005, 2014), bio.const = FALSE,sel.years = c(2005, 2014), sel.const = FALSE,Fcv=Fcv, Fphi=Fphi,Btrigger = 0,Blim=Blim,Bpa=Bpa,Fscan = seq(0,1.2,len=61),verbose=FALSE)

png("R_SSB_saithe_full time period2015.png")
eqsr_plot(FIT,n=2e4)
dev.off()

png("plot_saithe_full time period no Btrigger but Fcv and Fphi.png")
eqsim_plot(SIM1,catch="FALSE")
dev.off()

Coby.fit(SIM1,outfile='saithe full time series no Btrigger but Fcv and Fphi')

## using SIM101: find p05 column, run down to SSB matching FMSY - get biomasss
## FMSY=0.35 (in the middle)
round(round(mean(SIM101$rbp[,4][SIM101$rbp$variable=='Spawning stock biomass'&SIM101$rbp$Ftarget ==0.34],SIM101$rbp[,4][SIM101$rbp$variable=='Spawning stock biomass'&SIM101$rbp$Ftarget ==0.36])),-3)
##

######  with HCR --------------------------------------------------------------
## used to get Fp.05 and to get lower&upper range of FMSY with Btrigger
SIM2 <- eqsim_run(FIT,bio.years = c(2005, 2014), bio.const = FALSE,sel.years = c(2005, 2014), sel.const = FALSE, Fcv=Fcv, Fphi=Fphi,Btrigger = Btrigger,Blim=Blim,Bpa=Bpa,Fscan = seq(0,1.2,len=61),verbose=FALSE)

png("plot_saithe_full time period with Btrigger and Fcv and Fphi.png")
eqsim_plot(SIM2,catch="FALSE")
dev.off()

Coby.fit(SIM2,outfile='saithe full time series with Btrigger and Fcv and Fphi')

